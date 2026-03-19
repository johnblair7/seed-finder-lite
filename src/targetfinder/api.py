"""
High-level API for TargetFinder: run from Python/Jupyter, cache index, cluster-friendly.
"""

from __future__ import annotations

import csv
import os
import gzip
import hashlib
import shutil
import urllib.request
from pathlib import Path
from typing import List, Optional, Tuple, Union

from .search import check_bowtie_installed, run_bowtie_search
from .annotate import load_gtf_genes, annotate_off_targets
from .score import score_and_filter_off_targets
from .seed_finder import (
    find_seed_hits as _find_seed_hits,
    find_seed_hits_in_promoter_fasta as _find_seed_hits_promoter_fasta,
)
from .main import load_guides as _load_guides_from_path

# Default subdir for cached reference data (build once, reuse)
DEFAULT_REF_DIR = "refs"
DEFAULT_CACHE_DIR = Path.home() / ".cache" / "targetfinder"


def find_seed_hits(
    guide: str,
    genome_fa: Union[str, Path, None] = None,
    gtf_path: Union[str, Path, None] = None,
    min_seed_length: int = 8,
    tss_window_bp: int = 2000,
    ref_dir: Union[str, Path, None] = None,
    gtf_url: Union[str, None] = None,
    cache_dir: Union[str, Path, None] = None,
) -> List[dict]:
    """
    Direct genome scan for gRNA seed matches (no aligner).

    Returns all sites where the PAM-proximal seed (last n bases of the guide)
    perfectly matches the genome next to NGG, within ±tss_window_bp of a TSS.
    Each hit includes chr, start, end, strand, seed_length, pam, protospacer,
    guide_sequence, and ucsc_url.

    Provide either:
    - `ref_dir` (path to refs/ with genome.fa and gencode.gtf), or
    - `genome_fa` and `gtf_path`,
    - OR `gtf_url` (download+cache) to obtain `gtf_path`.

    Notes:
    - `genome_fa` can be the sentinel string `"ensembl"` to use Ensembl REST
      instead of a local FASTA (no local genome file needed).
    """

    cache_root = Path(cache_dir) if cache_dir is not None else DEFAULT_CACHE_DIR

    resolved_gtf_path: Optional[Path] = Path(gtf_path) if gtf_path is not None else None
    if resolved_gtf_path is None and ref_dir is not None:
        paths = get_ref_paths(Path(ref_dir))
        resolved_gtf_path = Path(paths["gtf"])

    if gtf_url is not None:
        url_str = str(gtf_url)
        h = hashlib.sha256(url_str.encode("utf-8")).hexdigest()[:16]
        gtf_cache_dir = cache_root / "gtf"
        gtf_cache_dir.mkdir(parents=True, exist_ok=True)

        is_gz = url_str.lower().endswith(".gz")
        cached_gz = gtf_cache_dir / f"{h}.gtf.gz"
        cached_gtf = gtf_cache_dir / f"{h}.gtf"

        # Download (or reuse) the gz file
        if is_gz:
            if not cached_gtf.exists():
                if not cached_gz.exists():
                    tmp = cached_gz.with_suffix(".tmp")
                    with urllib.request.urlopen(url_str, timeout=120) as resp:
                        with open(tmp, "wb") as f:
                            shutil.copyfileobj(resp, f)
                    tmp.replace(cached_gz)
                # Decompress to cached_gtf
                with gzip.open(cached_gz, "rb") as f_in, open(cached_gtf, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            resolved_gtf_path = cached_gtf
        else:
            # Plain gtf
            if not cached_gtf.exists():
                tmp = cached_gtf.with_suffix(".tmp")
                with urllib.request.urlopen(url_str, timeout=120) as resp:
                    with open(tmp, "wb") as f:
                        shutil.copyfileobj(resp, f)
                tmp.replace(cached_gtf)
            resolved_gtf_path = cached_gtf

    if genome_fa is None and ref_dir is not None:
        paths = get_ref_paths(Path(ref_dir))
        genome_fa = paths["genome_fa"]

    if not genome_fa or resolved_gtf_path is None:
        raise ValueError("Provide ref_dir, or genome_fa + gtf_path, or gtf_url (download+cache).")
    return _find_seed_hits(
        guide=guide,
        genome_fa=str(genome_fa),
        gtf_path=str(resolved_gtf_path),
        min_seed_length=min_seed_length,
        tss_window_bp=tss_window_bp,
    )


def find_seed_hits_promoter_fasta(
    guide: str,
    promoter_fa: Union[str, Path],
    gtf_path: Union[str, Path, None] = None,
    *,
    min_seed_length: int = 8,
    tss_window_bp: int = 2000,
    tss_table_path: Union[str, Path, None] = None,
) -> List[dict]:
    """
    Direct seed scan using a promoter-trimmed FASTA (headers like `prom|chr|s|e`).

    This bypasses creating TSS windows on the fly; it scans the contigs and
    remaps hits back to genomic coordinates, then filters to ±tss_window_bp
    of the closest TSS for annotation correctness.
    """
    promoter_fa_path = Path(promoter_fa)
    promoter_fa_to_use = promoter_fa_path

    # pyfaidx cannot read block-gz FASTA files reliably; decompress once into cache.
    if promoter_fa_path.suffix == ".gz":
        cache_dir = DEFAULT_CACHE_DIR / "promoter_fa"
        cache_dir.mkdir(parents=True, exist_ok=True)
        out_path = cache_dir / promoter_fa_path.stem  # keep the .fa suffix in the name
        if not out_path.exists():
            with gzip.open(promoter_fa_path, "rb") as f_in, open(out_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        promoter_fa_to_use = out_path

    return _find_seed_hits_promoter_fasta(
        guide=guide,
        promoter_fa=str(promoter_fa_to_use),
        gtf_path=str(gtf_path) if gtf_path is not None else None,
        tss_table_path=str(tss_table_path) if tss_table_path is not None else None,
        min_seed_length=min_seed_length,
        tss_window_bp=tss_window_bp,
    )


def get_ref_paths(
    ref_dir: Union[str, Path],
    *,
    index_name: str = "hg38",
    genome_name: str = "genome.fa",
    gtf_name: str = "gencode.gtf",
) -> dict:
    """
    Resolve paths for reference data under a single directory (cached index setup).

    Build the Bowtie2 index once in ref_dir; thereafter point all runs at ref_dir.
    """
    ref_dir = Path(ref_dir)
    return {
        "bowtie_index": str(ref_dir / index_name),
        "genome_fa": str(ref_dir / genome_name),
        "gtf": str(ref_dir / gtf_name),
    }


def load_guides(path: Union[str, Path]) -> List[Tuple[str, str]]:
    """Load (guide_ID, sequence) pairs from CSV/TSV. Public API."""
    from pathlib import Path
    return _load_guides_from_path(Path(path))


def run_off_target_analysis(
    guides: Union[str, Path, List[Tuple[str, str]]],
    *,
    bowtie_index: str,
    genome_fa: str,
    gtf: str,
    max_mismatches: int = 3,
    max_alignments_per_guide: int = 500,
    output_path: Union[str, Path, None] = None,
    return_dataframe: bool = True,
    exclude_same_locus: bool = True,
):
    """
    Run full off-target pipeline: align guides, filter by PAM, annotate with gene/TSS.

    Intended for use from a Jupyter notebook or scripts. The Bowtie2 index and
    genome/GTF paths can point to a cached refs directory so the index is built once.

    Parameters
    ----------
    guides : str, Path, or list of (guide_id, sequence)
        Input: path to CSV/TSV with guide_ID and sequence, or list of (id, seq) tuples.
    bowtie_index : str
        Path to Bowtie2 index (e.g. refs/hg38). Same index is reused across runs.
    genome_fa : str
        Genome FASTA used to build the index (for PAM extraction).
    gtf : str
        Gene annotations GTF (e.g. GENCODE).
    max_mismatches : int
        0–3.
    max_alignments_per_guide : int
        Bowtie2 -k; cap alignments per guide for speed. Use a large value or -a on cluster.
    output_path : str, Path, or None
        If set, write results CSV here.
    return_dataframe : bool
        If True and pandas is available, return a DataFrame; else list of dicts.
    exclude_same_locus : bool
        If True (default), drop hits in the same gene as the intended target (on-target).
        Set False to keep same-locus hits in the output with score 0.

    Returns
    -------
    pandas.DataFrame or list[dict]
        Off-target table with guide_id, chr, start, end, strand, mismatches,
        pam, protospacer, nearest_gene, tss_distance.
    """
    if isinstance(guides, (str, Path)):
        guide_list = load_guides(guides)
    else:
        guide_list = list(guides)

    if not guide_list:
        raise ValueError("No valid guides provided.")

    off_targets = run_bowtie_search(
        guide_list,
        bowtie_index,
        genome_fa,
        max_mismatches=max_mismatches,
        max_alignments_per_guide=max_alignments_per_guide,
    )
    if os.environ.get("TARGETFINDER_DEBUG"):
        import sys
        print(f"[TargetFinder] after Bowtie/parse: {len(off_targets)} records", file=sys.stderr)

    genes, tss_list = load_gtf_genes(gtf)
    results = annotate_off_targets(off_targets, genes, tss_list)
    if os.environ.get("TARGETFINDER_DEBUG"):
        import sys
        print(f"[TargetFinder] after annotate: {len(results)} records", file=sys.stderr)

    # Score and filter - ensure all records have required keys
    n_before_score = len(results)
    try:
        results = score_and_filter_off_targets(results, exclude_same_locus=exclude_same_locus)
        if os.environ.get("TARGETFINDER_DEBUG"):
            import sys
            print(f"[TargetFinder] after score/filter (exclude_same_locus={exclude_same_locus}): {len(results)} records", file=sys.stderr)
        # Verify scoring worked - check a sample record
        if results and len(results) > 0:
            sample = results[0]
            required_keys = ["score", "flag", "same_locus", "guide_sequence", "aligned_sequence"]
            missing = [k for k in required_keys if k not in sample]
            if missing:
                import warnings
                warnings.warn(f"Warning: Missing keys in scored results: {missing}. "
                            f"Sample record keys: {list(sample.keys())}")
    except Exception as e:
        import warnings
        warnings.warn(f"Error in scoring: {e}. Results may be missing score/flag columns.")
        # Ensure records have at least empty values for these keys
        for rec in results:
            rec.setdefault("score", "")
            rec.setdefault("flag", "")
            rec.setdefault("same_locus", "")
            rec.setdefault("guide_sequence", rec.get("guide_sequence", ""))
            rec.setdefault("aligned_sequence", rec.get("aligned_sequence", ""))

    if output_path is not None:
        _write_results_csv(results, Path(output_path))
    if os.environ.get("TARGETFINDER_DEBUG") and len(results) == 0 and n_before_score > 0:
        import sys
        print("[TargetFinder] All records were dropped by scoring/filter. Try exclude_same_locus=False to include on-target.", file=sys.stderr)

    if return_dataframe:
        try:
            import pandas as pd
            return pd.DataFrame(results)
        except ImportError:
            pass
    return results


def _write_results_csv(results: list, path: Path) -> None:
    cols = [
        "guide_id", "chr", "start", "end", "strand",
        "mismatches", "pam", "protospacer",
        "guide_sequence", "aligned_sequence",
        "nearest_gene", "tss_distance",
        "score", "flag", "same_locus", "bd_other_gene",
    ]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(results)


def build_bowtie_index(
    genome_fa: Union[str, Path],
    index_path: Union[str, Path],
) -> None:
    """
    Build Bowtie2 index once; reuse index_path for all future runs (cache).

    Run this once per genome (e.g. in refs/). Then use index_path in
    run_off_target_analysis(..., bowtie_index=index_path).
    """
    import subprocess
    genome_fa = Path(genome_fa)
    index_path = Path(index_path)
    if not genome_fa.exists():
        raise FileNotFoundError(f"Genome FASTA not found: {genome_fa}")
    subprocess.run(
        ["bowtie2-build", str(genome_fa), str(index_path)],
        check=True,
    )
