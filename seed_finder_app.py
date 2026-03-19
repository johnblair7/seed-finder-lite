"""
Seed Finder app — enter a 20 bp gRNA, get genomic coordinates and UCSC links
for all PAM-proximal seed matches in TSS regions (direct genome scan, no Bowtie).

Run from repo root:
  streamlit run seed_finder_app.py
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

import streamlit as st
from targetfinder.api import find_seed_hits, find_seed_hits_promoter_fasta, get_ref_paths


def main():
    st.set_page_config(page_title="Seed Finder", layout="wide")
    st.title("Seed Finder")
    st.markdown(
        "Enter a 20 bp gRNA sequence. The app scans the genome (TSS ± 2000 bp) for "
        "perfect **PAM-proximal seed** matches next to NGG and returns coordinates + UCSC links."
    )

    use_ensembl = st.checkbox(
        "Use Ensembl genome online (no local FASTA)",
        value=False,
        help="Fetch sequence from Ensembl REST API; only a local GTF is needed (for TSS regions).",
    )
    ref_dir = st.text_input(
        "Reference directory (for GTF; also for genome.fa if not using Ensembl)",
        value=str(ROOT / "refs"),
        help="Folder containing gencode.gtf (and genome.fa unless using Ensembl).",
    )
    ref_path = Path(ref_dir)
    paths = get_ref_paths(ref_path) if ref_path.is_dir() else {}
    if not ref_path.is_dir():
        st.warning(f"Directory not found: {ref_path}")
    elif paths:
        if not use_ensembl and not Path(paths["genome_fa"]).exists():
            st.warning(f"Genome not found: {paths['genome_fa']}")
        if not Path(paths["gtf"]).exists():
            st.warning(f"GTF not found: {paths['gtf']}")

    guide = st.text_input(
        "gRNA sequence (20 bp)",
        value="",
        placeholder="e.g. CCGGCGCCCGCAGAGCCCCG",
        max_chars=20,
    ).strip().upper()

    col1, col2 = st.columns(2)
    with col1:
        min_seed = st.slider("Min seed length", min_value=5, max_value=20, value=8)
    with col2:
        tss_window = st.number_input("TSS window (bp)", min_value=500, max_value=10000, value=2000)

    use_remote_gtf = st.checkbox(
        "Use remote GTF URL (download+cache)",
        value=False,
        help="If enabled, the app downloads the GTF once and reuses it from cache.",
    )
    default_ensembl_gtf = "https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz"
    gtf_url = st.text_input(
        "GTF URL",
        value=default_ensembl_gtf,
        help="Used only when 'Use remote GTF URL' is enabled.",
    )

    default_promoter_fa = str(ROOT / "refs" / "promoter_trim_2000_pad23.fa.gz")
    use_promoter_fa_mode = st.checkbox(
        "Promoter FASTA mode (scan trimmed promoter FASTA)",
        value=True,
        help="Scans a pre-trimmed FASTA containing only ±TSS windows (faster). Requires a promoter FASTA + GTF.",
    )
    promoter_fa_path = st.text_input(
        "Promoter FASTA path",
        value=default_promoter_fa,
        help="Headers must be like: prom|chr16|start|end (this app ships one in refs/). .fa.gz is OK (decompressed to cache).",
        disabled=not use_promoter_fa_mode,
    )

    default_tss_table = str(ROOT / "refs" / "tss_table.tsv")
    tss_table_path = st.text_input(
        "TSS table path (lite)",
        value=default_tss_table,
        help="TSV with columns: chr, tss, gene_name. Used in promoter FASTA mode to avoid shipping the full GTF.",
        disabled=not use_promoter_fa_mode,
    )

    if st.button("Find seed matches"):
        if not guide or len(guide) != 20:
            st.error("Please enter a 20 bp sequence.")
            return
        if not all(c in "ACGT" for c in guide):
            st.error("Sequence must be ACGT only.")
            return
        if use_promoter_fa_mode:
            if not promoter_fa_path.strip():
                st.error("Provide a promoter FASTA path.")
                return
            if not Path(promoter_fa_path).exists():
                st.error(f"Promoter FASTA not found: {promoter_fa_path}")
                return
            # Genome FASTA not required in this mode.
        elif not use_ensembl:
            # Local genome required
            if not ref_path.is_dir() or not paths:
                st.error("Reference directory not found. Add genome.fa to refs/ (or enable Ensembl genome online).")
                return
            if not Path(paths.get("genome_fa", "")).exists():
                st.error("Genome FASTA not found in refs/. Use 'Use Ensembl genome online' or add genome.fa to refs.")
                return
        if not use_promoter_fa_mode and not use_remote_gtf:
            # Local GTF required
            if not ref_path.is_dir() or not paths or not Path(paths.get("gtf", "")).exists():
                st.error("GTF not found. Either add gencode.gtf to refs/ or enable 'Use remote GTF URL'.")
                return
        if (not use_promoter_fa_mode) and use_remote_gtf and not gtf_url.strip():
            st.error("Provide a non-empty GTF URL.")
            return

        with st.spinner(
            "Scanning TSS regions..."
            + (" (Ensembl)" if use_ensembl else "")
            + (" (promoter FASTA mode)" if use_promoter_fa_mode else "")
        ):
            try:
                if use_promoter_fa_mode:
                    if not tss_table_path.strip() or not Path(tss_table_path).exists():
                        st.error(f"TSS table not found: {tss_table_path}")
                        return
                    results = find_seed_hits_promoter_fasta(
                        guide,
                        promoter_fa=promoter_fa_path,
                        gtf_path=None,
                        tss_table_path=tss_table_path,
                        min_seed_length=min_seed,
                        tss_window_bp=tss_window,
                    )
                else:
                    results = find_seed_hits(
                        guide,
                        genome_fa="ensembl" if use_ensembl else None,
                        gtf_path=(paths["gtf"] if (not use_remote_gtf and paths) else None),
                        gtf_url=(gtf_url if use_remote_gtf else None),
                        ref_dir=(ref_path if not use_ensembl else None),
                        min_seed_length=min_seed,
                        tss_window_bp=tss_window,
                    )
            except Exception as e:
                st.exception(e)
                return

        st.success(f"Found **{len(results)}** seed match(es).")

        if not results:
            st.info("No matches in TSS regions. Try a larger TSS window or a different guide.")
            return

        st.subheader("UCSC links")
        for i, r in enumerate(results, 1):
            url = r.get("ucsc_url", "")
            chrom = r.get("chr", "")
            start = r.get("start", "")
            end = r.get("end", "")
            strand = r.get("strand", "")
            seed_len = r.get("seed_length", "")
            gene = r.get("nearest_tss_gene", "")
            dist = r.get("tss_distance", "")
            st.markdown(
                f"{i}. **protospacer {chrom}:{start}-{end}** ({strand}, seed={seed_len}, {gene} ±{dist} bp) — [Open in UCSC]({url})"
            )

        st.divider()
        st.subheader("Table")
        import pandas as pd
        df = pd.DataFrame(results)
        st.dataframe(df, use_container_width=True, hide_index=True)


if __name__ == "__main__":
    main()
