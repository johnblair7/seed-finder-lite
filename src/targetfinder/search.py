"""Off-target search using Bowtie2."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from .seq_utils import reverse_complement


def check_bowtie_installed() -> bool:
    """Return True if Bowtie2 is available."""
    return shutil.which("bowtie2") is not None


def run_bowtie_search(
    guides: list[tuple[str, str]],
    bowtie_index: str,
    genome_fasta: str,
    max_mismatches: int = 3,
    max_alignments_per_guide: int = 500,
) -> list[dict]:
    """
    Search for off-targets using Bowtie2.

    Args:
        guides: List of (guide_id, sequence) tuples.
        bowtie_index: Path to Bowtie2 index (e.g. /path/to/hg38).
        genome_fasta: Path to genome FASTA (for PAM extraction).
        max_mismatches: Max mismatches (0-3).
        max_alignments_per_guide: Bowtie2 -k; max alignments per guide (use large value for cluster).

    Returns:
        List of off-target records with chr, start, end, strand, mismatches,
        guide_id, pam, protospacer.
    """
    try:
        import pyfaidx
    except ImportError:
        raise ImportError("pyfaidx is required for genome access. pip install pyfaidx")

    genome = pyfaidx.Fasta(genome_fasta)
    results = []

    with tempfile.TemporaryDirectory() as tmpdir:
        tmppath = Path(tmpdir)
        fasta_path = tmppath / "guides.fa"
        sam_path = tmppath / "alignments.sam"

        with open(fasta_path, "w") as f:
            for gid, seq in guides:
                seq = seq.upper().strip()
                if len(seq) != 20:
                    continue
                protospacer = reverse_complement(seq)
                f.write(f">{gid}\n{protospacer}\n")

        cmd = [
            "bowtie2",
            "-f",
            "-x", bowtie_index,
            "-U", str(fasta_path),
            "-S", str(sam_path),
            "-k", str(max_alignments_per_guide),
            "-N", "1",  # 1 mismatch in seed
            "-L", "10",  # seed length
            "--score-min", "L,0,-0.4",  # permissive for up to ~3 mismatches
            "--mp", "1,1",
        ]
        subprocess.run(cmd, check=True, capture_output=True, text=True)

        # Build guide_id -> original sequence mapping
        guide_seq_map = {gid: seq for gid, seq in guides if len(seq.strip().upper()) == 20}
        
        for rec in _parse_sam(str(sam_path), genome, max_mismatches, guide_seq_map):
            results.append(rec)

    genome.close()
    return results


def _parse_sam(sam_path: str, genome, max_mismatches: int = 3, guide_seq_map: dict[str, str] | None = None) -> list[dict]:
    """Parse Bowtie SAM output and filter for PAM-adjacent hits."""
    out = []
    debug = os.environ.get("TARGETFINDER_DEBUG")
    n_sam, n_skip_ref, n_skip_pam, n_skip_nm = 0, 0, 0, 0
    with open(sam_path) as f:
        for line in f:
            if line.startswith("@"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 11:
                continue
            n_sam += 1
            qname = parts[0]
            flag = int(parts[1])
            ref = parts[2]
            pos_1based = int(parts[3])
            cigar = parts[5]
            seq = parts[9]

            if ref == "*" or ref.startswith("chrUn") or "alt" in ref.lower():
                n_skip_ref += 1
                continue

            strand = "-" if (flag & 16) else "+"
            start_0 = pos_1based - 1
            cigar_parsed = _parse_cigar(cigar, len(seq))
            if cigar_parsed is None:
                ref_span = len(seq)
                indel_read_positions: list[tuple[int, str]] = []
            else:
                ref_span, read_span, indel_read_positions = cigar_parsed
                if read_span != len(seq):
                    ref_span = len(seq)
                    indel_read_positions = []
            end_0 = start_0 + ref_span

            pam, pam_fail_reason = _get_pam(genome, ref, start_0, end_0, strand)
            if pam is None:
                n_skip_pam += 1
                if debug:
                    print(f"[TargetFinder] skip PAM ref={ref} pos={pos_1based} guide={qname}: {pam_fail_reason}", file=sys.stderr)
                continue

            nm = _get_nm(parts)
            if nm is None:
                nm = _count_mismatches_from_md(parts)
            if nm > max_mismatches:
                n_skip_nm += 1
                continue

            # Mismatch positions 1-based from PAM (1 = PAM-proximal) for scoring.
            mismatch_positions_from_pam = []
            for p in parts[11:]:
                if p.startswith("MD:Z:"):
                    md_str = p[5:]
                    read_mis = _parse_md_mismatch_positions(len(seq), md_str)
                    mismatch_positions_from_pam = _mismatch_positions_from_pam(
                        len(seq), strand, read_mis
                    )
                    break
            # Add indel positions (single-base insertions/shifts) so they are scored like mismatches
            read_len = len(seq)
            for read_pos, _kind in indel_read_positions:
                if 0 <= read_pos < read_len:
                    from_pam = read_len - read_pos
                    if from_pam not in mismatch_positions_from_pam:
                        mismatch_positions_from_pam.append(from_pam)
            mismatch_positions_from_pam = sorted(mismatch_positions_from_pam)

            # Length of uninterrupted perfect match directly adjacent to PAM (position 1 = PAM-proximal)
            seed_length = 20
            for pos_from_pam in mismatch_positions_from_pam:
                if 1 <= pos_from_pam <= 20:
                    seed_length = pos_from_pam - 1
                    break

            # Extract reference sequence (genome sequence) at this site (full ref span for indels)
            aligned_sequence = _get_reference_sequence(genome, ref, start_0, end_0, strand)
            
            # Get original guide sequence
            guide_sequence = ""
            if guide_seq_map:
                guide_sequence = guide_seq_map.get(qname, "")

            out.append({
                "guide_id": qname,
                "chr": ref,
                "start": start_0 + 1,
                "end": end_0,
                "strand": strand,
                "mismatches": nm,
                "pam": pam,
                "protospacer": seq,
                "guide_sequence": guide_sequence,
                "aligned_sequence": aligned_sequence,
                "mismatch_positions_from_pam": mismatch_positions_from_pam,
                "seed_length": seed_length,
            })
    if debug:
        print(f"[TargetFinder] SAM: {n_sam} alignments, skip_ref={n_skip_ref} skip_pam={n_skip_pam} skip_nm={n_skip_nm} -> kept {len(out)}", file=sys.stderr)
    return out


def _chrom_keys(chrom: str):
    """Yield chrom and alternate (SAM may use '9', FASTA may use 'chr9')."""
    yield chrom
    if chrom.startswith("chr"):
        yield chrom.lstrip("chr") or chrom
    else:
        yield "chr" + chrom


def _get_reference_sequence(genome, chrom: str, start: int, end: int, strand: str) -> str:
    """
    Extract the reference (genome) sequence at the aligned site, in guide orientation.
    """
    for ch in _chrom_keys(chrom):
        try:
            ref_seq = str(genome[ch][start:end]).upper()
            if strand == "-":
                ref_seq = reverse_complement(ref_seq)
            return ref_seq
        except (IndexError, KeyError, TypeError):
            continue
    return ""


def _get_pam(genome, chrom: str, start: int, end: int, strand: str) -> tuple[str | None, str | None]:
    """Extract 3bp PAM (SpCas9) adjacent to the protospacer. Returns (pam, None) or (None, reason).
    Accepted: NGG only.
    """
    last_error = "chrom not in genome"
    for ch in _chrom_keys(chrom):
        try:
            # PAM is the 3 bp 3' of the protospacer on the *target* strand. We send RC(guide) to
            # Bowtie, so: SAM "+" = read on plus → target on minus → PAM at [start-3, start), ref
            # has plus-strand there (CCN), RC to get NGG. SAM "-" = read on minus → target on plus
            # → PAM at [end, end+3], read as-is (NGG).
            if strand == "+":
                raw = str(genome[ch][start - 3 : start]).upper()
                pam = reverse_complement(raw)
            else:
                pam = str(genome[ch][end : end + 3]).upper()
            if len(pam) != 3:
                return (None, f"PAM slice len={len(pam)} (chrom={ch})")
            # SpCas9: NGG only
            if pam[2] == "G" and pam[1] == "G":
                return (pam, None)
            return (None, f"PAM not NGG: '{pam}' (chrom={ch})")
        except (KeyError, TypeError) as e:
            last_error = f"chrom '{ch}' not in genome: {e}"
            continue
        except IndexError as e:
            last_error = f"index error for chrom={ch} start={start} end={end}: {e}"
            continue
    return (None, last_error)


def _get_nm(parts: list[str]) -> int | None:
    """Extract NM or XM tag from SAM optional fields."""
    for p in parts[11:]:
        if p.startswith("NM:i:"):
            return int(p.split(":")[2])
        if p.startswith("XM:i:"):
            return int(p.split(":")[2])
    return None


def _parse_md_mismatch_positions(read_len: int, md_str: str) -> list[int]:
    """
    Parse SAM MD:Z tag to get 0-based read positions of mismatches.
    MD format: e.g. "10A5G2" = 10 match, ref A (mismatch), 5 match, ref G (mismatch), 2 match.
    """
    positions: list[int] = []
    read_pos = 0
    i = 0
    while i < len(md_str) and read_pos < read_len:
        if md_str[i].isdigit():
            j = i
            while j < len(md_str) and md_str[j].isdigit():
                j += 1
            run = int(md_str[i:j])
            read_pos += run
            i = j
        elif md_str[i] in "ACGTN":
            positions.append(read_pos)
            read_pos += 1
            i += 1
        else:
            i += 1
    return positions


def _mismatch_positions_from_pam(
    read_len: int, strand: str, read_0based_positions: list[int]
) -> list[int]:
    """
    Convert 0-based read positions to 1-based positions from PAM.
    Position 1 = PAM-proximal (base adjacent to PAM), 20 = PAM-distal.
    The read's 3' end is always PAM-adjacent (both strands); use read_len - p so
    that position is from_pam 1 (fixes opposite-strand scoring).
    """
    out = []
    for p in read_0based_positions:
        if p < 0 or p >= read_len:
            continue
        # PAM is 3' of the protospacer; read 3' = last base of read → from_pam 1 for both strands
        from_pam = read_len - p
        out.append(from_pam)
    return sorted(out)


def _parse_cigar(cigar: str, seq_len: int) -> tuple[int, int, list[tuple[int, str]]] | None:
    """
    Parse CIGAR string. Returns (ref_span, read_span, indel_positions) or None if invalid.
    indel_positions: list of (read_pos_0based, 'I'|'D'). For D, read_pos is the position
    of the read base immediately before the deletion (so the gap is "after" that base).
    """
    if not cigar or cigar == "*":
        return None
    ref_span = 0
    read_span = 0
    indels: list[tuple[int, str]] = []
    i = 0
    while i < len(cigar):
        j = i
        while j < len(cigar) and cigar[j].isdigit():
            j += 1
        if j == i:
            return None
        length = int(cigar[i:j])
        if j >= len(cigar):
            return None
        op = cigar[j]
        i = j + 1
        if op == "M":
            ref_span += length
            read_span += length
        elif op == "I":
            for k in range(length):
                indels.append((read_span + k, "I"))
            read_span += length
        elif op == "D":
            for _ in range(length):
                indels.append((read_span, "D"))  # D is after current read position
            ref_span += length
        elif op in "SHP":
            if op == "S":
                read_span += length
            # H and P don't consume read or ref
        elif op == "N":
            ref_span += length
        else:
            return None
    return (ref_span, read_span, indels)


def _count_mismatches_from_md(parts: list[str]) -> int:
    """Fallback: approximate mismatches if NM not present."""
    for p in parts[11:]:
        if p.startswith("MD:Z:"):
            md = p[5:]
            return sum(1 for c in md if c in "ACGT")
    return 0
