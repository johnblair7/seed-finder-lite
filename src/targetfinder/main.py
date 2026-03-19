"""Main entry point for TargetFinder - CRISPR off-target analysis."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def main() -> None:
    """Run TargetFinder off-target analysis or seed finder."""
    if len(sys.argv) >= 2 and sys.argv[1] == "seed-find":
        _main_seed_find(sys.argv[2:])
        return

    parser = argparse.ArgumentParser(
        description="Find and annotate CRISPR guide RNA off-targets in the human genome."
    )
    parser.add_argument(
        "input",
        type=Path,
        help="Input CSV/TSV with columns guide_ID and sequence (20 bp each)",
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=None,
        help="Output CSV (default: stdout)",
    )
    parser.add_argument(
        "--bowtie-index",
        type=Path,
        required=True,
        help="Path to Bowtie index (e.g. refs/hg38)",
    )
    parser.add_argument(
        "--genome",
        type=Path,
        required=True,
        help="Path to genome FASTA (for PAM extraction)",
    )
    parser.add_argument(
        "--gtf",
        type=Path,
        required=True,
        help="Path to gene annotations GTF (e.g. GENCODE)",
    )
    parser.add_argument(
        "-n", "--max-mismatches",
        type=int,
        default=3,
        choices=[0, 1, 2, 3],
        help="Max mismatches (default: 3)",
    )
    parser.add_argument(
        "-k", "--max-alignments",
        type=int,
        default=500,
        help="Max alignments per guide (default: 500; use large value on cluster)",
    )
    parser.add_argument(
        "--include-same-locus",
        action="store_true",
        help="Include on-target/same-gene hits in output (default: exclude them)",
    )
    args = parser.parse_args()

    from .search import check_bowtie_installed, run_bowtie_search
    from .annotate import load_gtf_genes, annotate_off_targets
    from .score import score_and_filter_off_targets

    if not check_bowtie_installed():
        print("Error: Bowtie2 is not installed or not in PATH.", file=sys.stderr)
        print("Install with: brew install bowtie2", file=sys.stderr)
        sys.exit(1)

    guides = load_guides(args.input)
    if not guides:
        print("Error: No valid guides found in input.", file=sys.stderr)
        sys.exit(1)
    print(f"Loaded {len(guides)} guides.", file=sys.stderr)

    print("Searching for off-targets...", file=sys.stderr)
    off_targets = run_bowtie_search(
        guides,
        str(args.bowtie_index),
        str(args.genome),
        max_mismatches=args.max_mismatches,
        max_alignments_per_guide=args.max_alignments,
    )
    print(f"Found {len(off_targets)} off-target sites.", file=sys.stderr)

    print("Loading gene annotations...", file=sys.stderr)
    genes, tss_list = load_gtf_genes(str(args.gtf))
    print(f"Loaded {len(genes)} genes.", file=sys.stderr)

    print("Annotating with nearest gene and TSS distance...", file=sys.stderr)
    results = annotate_off_targets(off_targets, genes, tss_list)
    print("Scoring and filtering off-targets...", file=sys.stderr)
    try:
        results = score_and_filter_off_targets(
            results, exclude_same_locus=not args.include_same_locus
        )
        # Verify scoring worked
        if results and len(results) > 0:
            sample = results[0]
            required_keys = ["score", "flag", "same_locus", "guide_sequence", "aligned_sequence"]
            missing = [k for k in required_keys if k not in sample]
            if missing:
                print(f"Warning: Missing keys in scored results: {missing}", file=sys.stderr)
                print(f"Sample record keys: {list(sample.keys())}", file=sys.stderr)
        print(f"Reporting {len(results)} off-target sites (after seed filter).", file=sys.stderr)
    except ImportError as e:
        print(f"ERROR: Cannot import scoring module: {e}", file=sys.stderr)
        print("Make sure you ran: pip install -e .", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"ERROR in scoring: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc(file=sys.stderr)
        # Ensure records have at least empty values
        for rec in results:
            rec.setdefault("score", "")
            rec.setdefault("flag", "")
            rec.setdefault("same_locus", "")
            rec.setdefault("guide_sequence", "")
            rec.setdefault("aligned_sequence", "")
        print(f"Continuing with {len(results)} sites (scoring failed).", file=sys.stderr)

    _write_output(results, args.output)


def _main_seed_find(argv: list[str]) -> None:
    """CLI for direct seed-finder (no aligner)."""
    parser = argparse.ArgumentParser(
        description="Find gRNA seed matches in TSS regions (direct genome scan, no Bowtie)."
    )
    parser.add_argument(
        "guide",
        type=str,
        help="20 bp gRNA sequence (ACGT)",
    )
    parser.add_argument(
        "--refs",
        type=Path,
        default=Path("refs"),
        help="Reference directory with genome.fa and GTF (default: refs)",
    )
    parser.add_argument(
        "--min-seed-length",
        type=int,
        default=8,
        metavar="N",
        help="Minimum PAM-proximal seed length (default: 8)",
    )
    parser.add_argument(
        "--tss-window",
        type=int,
        default=2000,
        metavar="BP",
        help="Max distance from TSS in bp (default: 2000)",
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=None,
        help="Output CSV (default: print to stdout)",
    )
    args = parser.parse_args(argv)

    from .api import find_seed_hits, get_ref_paths

    paths = get_ref_paths(args.refs)
    genome_fa = paths["genome_fa"]
    gtf_path = paths["gtf"]
    if not Path(genome_fa).exists():
        print(f"Error: Genome not found: {genome_fa}", file=sys.stderr)
        sys.exit(1)
    if not Path(gtf_path).exists():
        print(f"Error: GTF not found: {gtf_path}", file=sys.stderr)
        sys.exit(1)

    print("Scanning TSS regions for seed matches...", file=sys.stderr)
    results = find_seed_hits(
        args.guide.strip().upper(),
        genome_fa=genome_fa,
        gtf_path=gtf_path,
        min_seed_length=args.min_seed_length,
        tss_window_bp=args.tss_window,
    )
    print(f"Found {len(results)} seed match(es).", file=sys.stderr)

    cols = ["chr", "start", "end", "strand", "seed_length", "nearest_tss_gene", "tss_distance", "pam", "protospacer", "guide_sequence", "ucsc_url"]
    out = open(args.output, "w", newline="") if args.output else sys.stdout
    writer = csv.DictWriter(out, fieldnames=cols, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(results)
    if args.output:
        out.close()


def load_guides(path: Path) -> list[tuple[str, str]]:
    """Load guide_ID and sequence from CSV/TSV."""
    guides = []
    path = Path(path)
    delim = "\t" if path.suffix in (".tsv", ".txt") else ","
    with open(path) as f:
        reader = csv.DictReader(f, delimiter=delim)
        cols = reader.fieldnames or []
        id_col = _find_col(cols, ["guide_id", "guide_ID", "id", "ID"])
        seq_col = _find_col(cols, ["sequence", "seq", "spacer"])
        if not id_col or not seq_col:
            raise ValueError(
                f"Input must have guide_ID and sequence columns. Found: {cols}"
            )
        for row in reader:
            gid = (row.get(id_col) or "").strip()
            seq = (row.get(seq_col) or "").strip().upper()
            if gid and seq and len(seq) == 20 and all(c in "ACGT" for c in seq):
                guides.append((gid, seq))
    return guides


def _find_col(columns: list[str], candidates: list[str]) -> str | None:
    for c in candidates:
        for col in columns:
            clean = col.strip().lower().lstrip("\ufeff")
            if clean == c.lower():
                return col
    return None


def _write_output(results: list[dict], path: Path | None) -> None:
    cols = [
        "guide_id", "chr", "start", "end", "strand",
        "mismatches", "pam", "protospacer",
        "guide_sequence", "aligned_sequence",
        "nearest_gene", "tss_distance",
        "score", "flag", "same_locus", "bd_other_gene",
    ]
    out = open(path, "w", newline="") if path else sys.stdout
    writer = csv.DictWriter(out, fieldnames=cols, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(results)
    if path:
        out.close()


if __name__ == "__main__":
    main()
