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
from targetfinder.api import find_seed_hits_promoter_fasta


def main():
    st.set_page_config(page_title="Seed Finder", layout="wide")
    st.title("Seed Finder")
    st.markdown(
        "Enter a 20 bp gRNA sequence. This **lite** version scans pre-trimmed TSS windows "
        "from the shipped promoter FASTA and returns PAM-proximal seed match coordinates + UCSC links."
    )

    # Lite mode is promoter FASTA-only (no Ensembl, no remote GTF).
    promoter_fa_path = ROOT / "refs" / "promoter_trim_2000_pad23.fa.gz"
    tss_table_path = ROOT / "refs" / "tss_table.tsv"
    if not promoter_fa_path.exists():
        st.error(f"Missing promoter FASTA: {promoter_fa_path}")
        return
    if not tss_table_path.exists():
        st.error(f"Missing TSS table: {tss_table_path}")
        return

    guide = st.text_input(
        "gRNA sequence (20 bp)",
        value="",
        placeholder="e.g. CCGGCGCCCGCAGAGCCCCG",
        max_chars=20,
    ).strip().upper()

    col1, col2 = st.columns(2)
    with col1:
        min_seed = st.slider("Min seed length", min_value=5, max_value=20, value=8)

    if st.button("Find seed matches"):
        if not guide or len(guide) != 20:
            st.error("Please enter a 20 bp sequence.")
            return
        if not all(c in "ACGT" for c in guide):
            st.error("Sequence must be ACGT only.")
            return

        # Promoter FASTA mode: only the shipped assets are used.

        with st.spinner(
            "Scanning TSS regions (lite promoter FASTA)..."
        ):
            try:
                results = find_seed_hits_promoter_fasta(
                    guide,
                    promoter_fa=str(promoter_fa_path),
                    gtf_path=None,
                    tss_table_path=str(tss_table_path),
                    min_seed_length=min_seed,
                    tss_window_bp=1000,
                )
            except Exception as e:
                st.exception(e)
                return

        st.success(f"Found **{len(results)}** seed match(es).")

        if not results:
            st.info("No matches in the default TSS window (1000 bp) for this guide.")
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
