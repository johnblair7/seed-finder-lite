"""
Direct genome scan for gRNA seed matches (no aligner).

Finds all sites where the PAM-proximal seed (last n bases of the guide) perfectly
matches the genome, adjacent to NGG (plus strand) or CCN (minus strand), restricted
to ±tss_window_bp of transcriptional start sites.
"""

from __future__ import annotations

import time
from functools import lru_cache
from pathlib import Path
from typing import List, Optional, Tuple

from .seq_utils import reverse_complement

# Use Ensembl REST API when genome_fa is this sentinel or an Ensembl base URL
ENSEMBL_GENOME = "ensembl"
ENSEMBL_API_BASE = "https://rest.ensembl.org"


def _chrom_to_ensembl(chrom: str) -> str:
    """Ensembl uses '1', 'X' etc.; strip 'chr' prefix."""
    return chrom[3:] if chrom.startswith("chr") else chrom


def _fetch_ensembl_sequence(
    chrom: str,
    start_1: int,
    end_1: int,
    base_url: str = ENSEMBL_API_BASE,
) -> Optional[str]:
    """
    Fetch genomic sequence from Ensembl REST API (1-based inclusive).
    Returns None on failure. Max 10 Mb per request.
    """
    try:
        import urllib.request
    except ImportError:
        return None
    region = f"{_chrom_to_ensembl(chrom)}:{start_1}..{end_1}"
    url = f"{base_url.rstrip('/')}/sequence/region/human/{region}?content-type=text/plain"
    try:
        req = urllib.request.Request(url, headers={"Content-Type": "text/plain"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            raw = resp.read().decode("utf-8")
    except Exception:
        return None
    # Response is FASTA: first line ">...", rest is sequence possibly with newlines
    lines = raw.strip().splitlines()
    seq = "".join(line.strip().upper() for line in lines if line and not line.startswith(">"))
    return seq if seq and all(c in "ACGTN" for c in seq) else None


def _load_tss_with_genes(gtf_path: str) -> List[Tuple[str, int, str]]:
    """Load (chrom, tss_1based, gene_name) from GTF gene entries."""
    out: List[Tuple[str, int, str]] = []
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            chrom = parts[0]
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            tss = start if strand == "+" else end
            attr = parts[8]
            gene_name = ""
            for item in attr.split(";"):
                item = item.strip()
                if item.lower().startswith("gene_name "):
                    gene_name = item.split(" ", 1)[1].strip().strip('"')
                    break
            out.append((chrom, tss, gene_name or "?"))
    return out


@lru_cache(maxsize=4)
def _load_tss_table(tsv_path: str) -> List[Tuple[str, int, str]]:
    """
    Load a precomputed TSS table (lite mode).

    Expected TSV columns:
      chr<TAB>tss<TAB>gene_name
    Optional header line is ignored.
    """
    out: List[Tuple[str, int, str]] = []
    with open(tsv_path) as f:
        first = f.readline()
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            tss = int(parts[1])
            gene_name = parts[2] or "?"
            out.append((chrom, tss, gene_name))
    return out


def _nearest_tss_gene(pos: int, chrom: str, tss_list: List[Tuple[str, int, str]]) -> Tuple[str, int]:
    """Return (gene_name, distance_bp) for the TSS closest to pos on chrom."""
    best_name = ""
    best_dist = float("inf")
    for c, tss, name in tss_list:
        if c != chrom:
            continue
        d = abs(pos - tss)
        if d < best_dist:
            best_dist = d
            best_name = name
    return (best_name or "", int(best_dist) if best_dist != float("inf") else 0)


def _chrom_keys(chrom: str):
    """Yield chrom and alternate (GTF may use 'chr1', FASTA may use '1')."""
    yield chrom
    if chrom.startswith("chr"):
        yield chrom[3:] or chrom
    else:
        yield "chr" + chrom


def load_tss_promoter_windows(
    gtf_path: str,
    tss_window_bp: int = 2000,
) -> List[Tuple[str, int, int]]:
    """
    Build merged promoter windows from GTF gene entries.

    Returns:
        List of (chrom, start_1based, end_1based) windows covering ±tss_window_bp
        around each TSS. Overlapping windows are merged per chromosome.
    """
    tss_list: List[Tuple[str, int]] = []  # (chrom, tss)
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            chrom = parts[0]
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            tss = start if strand == "+" else end
            tss_list.append((chrom, tss))

    # Group by chrom and merge overlapping intervals
    by_chr: dict[str, List[Tuple[int, int]]] = {}
    for chrom, tss in tss_list:
        s = max(1, tss - tss_window_bp)
        e = tss + tss_window_bp
        by_chr.setdefault(chrom, []).append((s, e))

    out: List[Tuple[str, int, int]] = []
    for chrom, intervals in by_chr.items():
        intervals.sort()
        cur_s, cur_e = intervals[0]
        for s, e in intervals[1:]:
            if s <= cur_e + 1:
                cur_e = max(cur_e, e)
            else:
                out.append((chrom, cur_s, cur_e))
                cur_s, cur_e = s, e
        out.append((chrom, cur_s, cur_e))
    return out


def _count_pam_proximal_matches(a: str, b: str) -> int:
    """Consecutive matches from PAM-proximal end (index 19 down). a and b length 20."""
    n = 0
    for i in range(19, -1, -1):
        if a[i] == b[i]:
            n += 1
        else:
            break
    return n


def _scan_sequence_for_hits(
    seq: str,
    chrom: str,
    s0: int,
    guide: str,
    guide_rc: str,
    guide_last_n: str,
    min_seed_length: int,
) -> List[dict]:
    """Scan one sequence window for plus/minus strand seed hits. s0 = 0-based genomic start of window."""
    hits: List[dict] = []
    if len(seq) < 23:
        return hits
    # PAM-proximal seed match (computed once; previously recomputed inside the loop)
    seed_rc = reverse_complement(guide_last_n)
    seed_start = 20 - min_seed_length
    for i in range(0, len(seq) - 23 + 1):
        # PAM NGG means positions i+21 and i+22 are both G
        if seq[i + 21] != "G" or seq[i + 22] != "G":
            continue
        # Perfect PAM-proximal seed match
        if seq[i + seed_start : i + 20] != seed_rc:
            continue
        prot = seq[i : i + 20]
        pam = seq[i + 20 : i + 23]
        seed_len = _count_pam_proximal_matches(prot, guide_rc)
        if seed_len < min_seed_length:
            continue
        start_1 = s0 + i + 1
        end_1 = start_1 + 19
        hits.append({
            "chr": chrom,
            "start": start_1,
            "end": end_1,
            # Coordinates are protospacer-only (PAM bases excluded) in genomic space.
            "protospacer_start": start_1,
            "protospacer_end": end_1,
            "pam_start": end_1 + 1,
            "pam_end": end_1 + 3,
            "strand": "-",
            "seed_length": seed_len,
            "pam": pam,
            "protospacer": prot,
            "guide_sequence": guide,
            "ucsc_url": _ucsc_url(chrom, start_1, end_1),
        })
    for i in range(0, len(seq) - 23 + 1):
        # Reverse-strand equivalent PAM is CCN: positions i and i+1 are both C
        if seq[i] != "C" or seq[i + 1] != "C":
            continue
        # Perfect seed match is on the 5' end of the protospacer for this orientation
        if seq[i + 3 : i + 3 + min_seed_length] != seed_rc:
            continue
        pam_ccn = seq[i : i + 3]
        prot_ref = seq[i + 3 : i + 23]
        prot_minus = reverse_complement(prot_ref)
        seed_len = _count_pam_proximal_matches(prot_minus, guide)
        if seed_len < min_seed_length:
            continue
        start_1 = s0 + i + 3 + 1
        end_1 = start_1 + 19
        hits.append({
            "chr": chrom,
            "start": start_1,
            "end": end_1,
            # For minus-target hits, PAM lies 3 bp upstream of the protospacer in genomic coordinates.
            "protospacer_start": start_1,
            "protospacer_end": end_1,
            "pam_start": start_1 - 3,
            "pam_end": start_1 - 1,
            "strand": "+",
            "seed_length": seed_len,
            "pam": reverse_complement(pam_ccn),
            "protospacer": reverse_complement(prot_ref),
            "guide_sequence": guide,
            "ucsc_url": _ucsc_url(chrom, start_1, end_1),
        })
    return hits


def find_seed_hits(
    guide: str,
    genome_fa: str,
    gtf_path: str,
    min_seed_length: int = 8,
    tss_window_bp: int = 2000,
) -> List[dict]:
    """
    Scan genome for perfect PAM-proximal seed matches within TSS promoter windows.

    Args:
        guide: 20 bp gRNA sequence (5' to 3').
        genome_fa: Path to reference FASTA, or "ensembl" to use Ensembl REST API (no local genome).
        gtf_path: Path to GTF (gene entries used for TSS).
        min_seed_length: Minimum number of PAM-proximal bases that must match (default 8).
        tss_window_bp: Max distance from TSS to include (default 2000).

    Returns:
        List of dicts with chr, start, end, strand, seed_length, pam, protospacer,
        guide_sequence, ucsc_url, nearest_tss_gene, and tss_distance (1-based coordinates).
    """
    guide = guide.strip().upper()
    if len(guide) != 20 or not all(c in "ACGT" for c in guide):
        raise ValueError("Guide must be 20 bp ACGT.")

    guide_rc = reverse_complement(guide)
    guide_last_n = guide[-min_seed_length:]
    windows = load_tss_promoter_windows(gtf_path, tss_window_bp)
    hits: List[dict] = []

    use_ensembl = (
        genome_fa.strip().lower() == ENSEMBL_GENOME
        or (isinstance(genome_fa, str) and "ensembl" in genome_fa.lower())
    )
    ensembl_base = genome_fa if (use_ensembl and "://" in genome_fa) else ENSEMBL_API_BASE

    if use_ensembl:
        for chrom, win_start, win_end in windows:
            seq = _fetch_ensembl_sequence(chrom, win_start, win_end + 3, ensembl_base)
            time.sleep(0.05)
            if not seq or len(seq) < 23:
                continue
            s0 = win_start - 1
            hits.extend(
                _scan_sequence_for_hits(
                    seq, chrom, s0, guide, guide_rc, guide_last_n, min_seed_length
                )
            )
    else:
        try:
            from pyfaidx import Fasta
        except ImportError:
            raise ImportError("pyfaidx is required for local genome. pip install pyfaidx")
        fa = Fasta(genome_fa)
        for chrom, win_start, win_end in windows:
            chrom_key = chrom
            for key in _chrom_keys(chrom):
                if key in fa:
                    chrom_key = key
                    break
            else:
                continue
            s0 = max(0, win_start - 1)
            e0 = min(win_end + 3, len(fa[chrom_key]))
            if e0 - s0 < 23:
                continue
            seq = str(fa[chrom_key][s0:e0]).upper()
            hits.extend(
                _scan_sequence_for_hits(
                    seq, chrom, s0, guide, guide_rc, guide_last_n, min_seed_length
                )
            )
        fa.close()

    # Annotate each hit with nearest TSS gene name and distance
    tss_with_genes = _load_tss_with_genes(gtf_path)
    for h in hits:
        pos = (h["start"] + h["end"]) // 2
        gene_name, dist = _nearest_tss_gene(pos, h["chr"], tss_with_genes)
        h["nearest_tss_gene"] = gene_name
        h["tss_distance"] = dist

    return hits


def _parse_promoter_contig_name(name: str) -> Tuple[str, int, int]:
    """
    Parse contig header from `promoter_trim_*` FASTA.

    Expected format: prom|{chrom}|{start_pad}|{end_pad}
    where start_pad/end_pad are 1-based inclusive genomic coordinates.
    """
    parts = (name or "").split("|")
    if len(parts) < 4 or parts[0] != "prom":
        raise ValueError(f"Unrecognized promoter contig name: {name!r}")
    chrom = parts[1]
    start_pad = int(parts[2])
    end_pad = int(parts[3])
    return chrom, start_pad, end_pad


def find_seed_hits_in_promoter_fasta(
    guide: str,
    promoter_fa: str,
    gtf_path: str = None,
    tss_table_path: str = None,
    min_seed_length: int = 8,
    tss_window_bp: int = 2000,
) -> List[dict]:
    """
    Scan a promoter-trimmed FASTA (e.g. extracted to TSS windows) for perfect
    PAM-proximal seed matches.

    The promoter FASTA must have headers like:
      prom|chr16|189257...|189267...
    so we can remap hits back to genomic coordinates.
    """
    try:
        from pyfaidx import Fasta
    except ImportError:
        raise ImportError("pyfaidx is required for seed finder. pip install pyfaidx")

    guide = guide.strip().upper()
    if len(guide) != 20 or not all(c in "ACGT" for c in guide):
        raise ValueError("Guide must be 20 bp ACGT.")

    guide_rc = reverse_complement(guide)
    guide_last_n = guide[-min_seed_length:]

    # Load TSS list once for distance + closest gene annotation.
    if tss_table_path is not None:
        tss_with_genes = _load_tss_table(tss_table_path)
    elif gtf_path is not None:
        tss_with_genes = _load_tss_with_genes(gtf_path)
    else:
        raise ValueError("Provide either gtf_path or tss_table_path.")

    fa = Fasta(promoter_fa)
    hits: List[dict] = []

    # pyfaidx exposes FASTA record names as keys.
    for contig in fa.keys():
        try:
            chrom, start_pad, end_pad = _parse_promoter_contig_name(contig)
        except Exception:
            continue

        # s0 is 0-based half-open base for `_scan_sequence_for_hits`
        s0 = start_pad - 1
        seq = str(fa[contig][:]).upper()
        if len(seq) < 23:
            continue

        local_hits = _scan_sequence_for_hits(
            seq=seq,
            chrom=chrom,
            s0=s0,
            guide=guide,
            guide_rc=guide_rc,
            guide_last_n=guide_last_n,
            min_seed_length=min_seed_length,
        )
        # Filter by true ±tss_window_bp requirement.
        for h in local_hits:
            pos = (h["start"] + h["end"]) // 2
            gene_name, dist = _nearest_tss_gene(pos, h["chr"], tss_with_genes)
            if dist > tss_window_bp:
                continue
            h["nearest_tss_gene"] = gene_name
            h["tss_distance"] = dist
            hits.append(h)

    fa.close()
    return hits


def _ucsc_url(chrom: str, start: int, end: int, db: str = "hg38") -> str:
    """UCSC Genome Browser link (1-based)."""
    c = chrom if str(chrom).startswith("chr") else f"chr{chrom}"
    return f"https://genome.ucsc.edu/cgi-bin/hgTracks?db={db}&position={c}:{start}-{end}"
