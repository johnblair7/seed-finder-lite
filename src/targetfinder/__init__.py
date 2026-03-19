"""TargetFinder - CRISPR guide off-target discovery and analysis."""

__version__ = "0.1.0"

from .api import (
    run_off_target_analysis,
    find_seed_hits,
    find_seed_hits_promoter_fasta,
    load_guides,
    get_ref_paths,
    build_bowtie_index,
)
from .cluster_utils import split_guides, merge_results

__all__ = [
    "run_off_target_analysis",
    "find_seed_hits",
    "find_seed_hits_promoter_fasta",
    "load_guides",
    "get_ref_paths",
    "build_bowtie_index",
    "split_guides",
    "merge_results",
]
