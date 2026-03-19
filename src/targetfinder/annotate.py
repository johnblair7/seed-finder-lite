"""Gene and TSS annotation for off-target sites."""

from __future__ import annotations

import os
import sys
from pathlib import Path


def load_gtf_genes(gtf_path: str) -> tuple[list[dict], list[dict]]:
    """
    Load genes and TSS positions from a GTF file.

    Returns:
        (genes, tss_list) where genes have chr, start, end, gene_id, gene_name, strand
        and tss_list has chr, tss, gene_id, gene_name, strand.
    """
    genes = []
    tss_list = []
    seen_genes = set()

    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != "gene":
                continue

            chrom = parts[0]
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attr = _parse_gtf_attrs(parts[8])

            gene_id = attr.get("gene_id", "")
            gene_name = attr.get("gene_name", gene_id)

            key = (chrom, start, end, gene_id)
            if key in seen_genes:
                continue
            seen_genes.add(key)

            genes.append({
                "chr": chrom,
                "start": start,
                "end": end,
                "gene_id": gene_id,
                "gene_name": gene_name,
                "strand": strand,
            })

            tss = start if strand == "+" else end
            tss_list.append({
                "chr": chrom,
                "tss": tss,
                "gene_id": gene_id,
                "gene_name": gene_name,
                "strand": strand,
            })

    return genes, tss_list


def _parse_gtf_attrs(s: str) -> dict[str, str]:
    out = {}
    for item in s.split(";"):
        item = item.strip()
        if not item:
            continue
        if " " in item:
            k, v = item.split(" ", 1)
            v = v.strip('"')
            out[k] = v
    return out


def annotate_off_targets(
    off_targets: list[dict],
    genes: list[dict],
    tss_list: list[dict],
) -> list[dict]:
    """
    Add nearest_gene and tss_distance to each off-target record.
    """
    # Index by chromosome
    genes_by_chr: dict[str, list[dict]] = {}
    for g in genes:
        c = g["chr"]
        if c not in genes_by_chr:
            genes_by_chr[c] = []
        genes_by_chr[c].append(g)

    tss_by_chr: dict[str, list[dict]] = {}
    for t in tss_list:
        c = t["chr"]
        if c not in tss_by_chr:
            tss_by_chr[c] = []
        tss_by_chr[c].append(t)

    out = []
    for rec in off_targets:
        chr_ = rec["chr"]
        chr_alt = "chr" + chr_ if not chr_.startswith("chr") else chr_.lstrip("chr")
        pos = (rec["start"] + rec["end"]) // 2

        genes_sub = genes_by_chr.get(chr_, []) or genes_by_chr.get(chr_alt, [])
        tss_sub = tss_by_chr.get(chr_, []) or tss_by_chr.get(chr_alt, [])

        nearest_gene, overlap = _nearest_gene(pos, genes_sub)
        # TSS distance is always for the gene the hit is in (nearest_gene). Never use the TARGET
        # gene's TSS for scoring—e.g. a hit in MTOR uses MTOR's TSS distance even if the guide
        # was designed for ANGPTL7 and ANGPTL7's TSS is closer.
        tss_dist = _tss_distance_for_gene(pos, tss_sub, nearest_gene) if nearest_gene else None
        if tss_dist is None:
            tss_dist = _nearest_tss_distance(pos, tss_sub)
        nearest_tss_gene, nearest_tss_d, nearest_tss_id, other_gene_within_2kb, min_dist_other, other_gene_name = _nearest_tss_and_other_gene_within_2kb(
            pos, tss_sub
        )

        rec = dict(rec)
        rec["nearest_gene"] = nearest_gene or ""
        rec["tss_distance"] = tss_dist if tss_dist is not None else ""
        rec["nearest_tss_gene"] = nearest_tss_gene or ""
        rec["nearest_tss_gene_id"] = nearest_tss_id or ""
        rec["other_gene_tss_within_2kb"] = other_gene_within_2kb
        rec["min_dist_other_gene_tss"] = min_dist_other if min_dist_other is not None else ""
        # BD other = the bidirectional partner, i.e. the gene with TSS within 2 kb that is NOT the target
        intended = (rec.get("guide_id") or "").split("_")[0].strip() if "_" in (rec.get("guide_id") or "") else ""
        intended_upper = intended.upper()
        if other_gene_within_2kb and intended_upper:
            partner = next(
                (g for g in (nearest_tss_gene, other_gene_name) if g and g.upper() != intended_upper),
                "",
            )
            rec["bd_other_gene"] = partner  # never the target; the true bidirectional partner
        else:
            rec["bd_other_gene"] = other_gene_name or ""
        if os.environ.get("TARGETFINDER_DEBUG") and rec.get("guide_id") == "SAXO1_0" and str(chr_) in ("9", "chr9"):
            print(
                f"[TargetFinder annotate] SAXO1_0 chr={chr_} pos={pos}: "
                f"nearest_gene={nearest_gene!r} tss_dist={tss_dist} "
                f"other_gene_tss_within_2kb={other_gene_within_2kb} bd_other_gene={other_gene_name!r} "
                f"#tss_sub={len(tss_sub)}",
                file=sys.stderr,
            )
        out.append(rec)

    return out


def _nearest_gene(pos: int, genes: list[dict]) -> tuple[str | None, bool]:
    """Find nearest gene. Return (gene_name, True if overlapping)."""
    if not genes:
        return None, False

    best_name = None
    best_dist = float("inf")
    overlapping = False

    for g in genes:
        if g["start"] <= pos <= g["end"]:
            overlapping = True
            return g["gene_name"], True
        d = min(abs(pos - g["start"]), abs(pos - g["end"]))
        if d < best_dist:
            best_dist = d
            best_name = g["gene_name"]

    return best_name, False


def _nearest_tss_distance(pos: int, tss_list: list[dict]) -> int | None:
    """Return distance to nearest TSS (bp)."""
    if not tss_list:
        return None
    return min(abs(pos - t["tss"]) for t in tss_list)


def _tss_distance_for_gene(pos: int, tss_list: list[dict], gene_name: str) -> int | None:
    """Return distance to the nearest TSS of the given gene (by gene_name). Used so
    scoring reflects the gene the hit is in (e.g. hit in MTOR uses MTOR's TSS distance),
    not a nested gene's TSS (e.g. ANGPTL7 TSS inside MTOR)."""
    if not gene_name or not tss_list:
        return None
    gene_upper = gene_name.upper()
    dists = [abs(pos - t["tss"]) for t in tss_list if (t.get("gene_name") or "").upper() == gene_upper]
    return min(dists) if dists else None


# Must match score.BD_MAX_TSS_DISTANCE_BP for BD flag
BD_TSS_DISTANCE_BP = 2000


def _nearest_tss_and_other_gene_within_2kb(
    pos: int, tss_list: list[dict]
) -> tuple[
    str | None, int | None, str | None,
    bool, int | None, str | None,
]:
    """
    Return (nearest_gene_name, nearest_dist, nearest_gene_id,
            other_gene_tss_within_2kb, min_dist_to_other_gene, other_gene_name).

    other_gene_name: name of the closest TSS (within 2 kb) that belongs to a
    *different* gene than the nearest; empty if none. "Different" = different
    gene_id, or different gene_name (catches bidirectional pairs that share
    gene_id in some GTFs, e.g. EXOSC8).
    """
    if not tss_list:
        return None, None, None, False, None, None
    with_dist = [
        (t["gene_name"], abs(pos - t["tss"]), t.get("gene_id", ""))
        for t in tss_list
    ]
    with_dist.sort(key=lambda x: x[1])
    nearest_name, nearest_d, nearest_id = with_dist[0]
    min_other_d = None
    other_gene_name = None
    for name, d, gid in with_dist[1:]:
        if d > BD_TSS_DISTANCE_BP:
            break
        # Other = different gene (by id or by name, so true BD pairs are caught)
        is_other = gid != nearest_id or (name and name != nearest_name)
        if is_other:
            if min_other_d is None or d < min_other_d:
                min_other_d = d
                other_gene_name = name
    other_gene_tss_within_2kb = min_other_d is not None
    return nearest_name, nearest_d, nearest_id, other_gene_tss_within_2kb, min_other_d, other_gene_name
