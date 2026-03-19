"""
Off-target scoring: TSS distance, PAM-proximal seed, same-locus, BD flag.
"""

from __future__ import annotations

import os
import sys

# Scoring weights (tunable)
TSS_POINTS_UNDER_200BP = 40
TSS_POINTS_AT_2000BP = 10
SEED_PERFECT_10_POINTS = 50
SEED_FIRST_MISMATCH_PENALTY_PER_POSITION = 5  # position 1 -> -5, position 10 -> -50
BD_BONUS = 100
TOTAL_MISMATCH_PENALTY = 5  # -5 per mismatch

# Hard filter: exclude if any mismatch in first 5 bases (PAM-proximal)
SEED_FILTER_FIRST_N_BASES = 5

# BD: on-target with another gene's TSS within this distance (bp)
BD_MAX_TSS_DISTANCE_BP = 2000


def _tss_score(tss_distance) -> float:
    """
    TSS distance score: closer = higher. Uses the hit's gene's TSS (set in annotate),
    never the target/intended gene's TSS.
    0-200 bp -> 40 pts; 200-2000 bp -> linear 10-40; >2000 -> 0.
    """
    if tss_distance is None or tss_distance == "":
        return 0.0
    try:
        d = int(tss_distance)
    except (TypeError, ValueError):
        return 0.0
    if d > 2000:
        return 0.0
    if d <= 200:
        return float(TSS_POINTS_UNDER_200BP)
    # Linear interpolation 200 -> 40, 2000 -> 10
    return TSS_POINTS_AT_2000BP + (TSS_POINTS_UNDER_200BP - TSS_POINTS_AT_2000BP) * (2000 - d) / (2000 - 200)


def _seed_score(mismatch_positions_from_pam: list[int]) -> tuple[float, bool]:
    """
    Seed (PAM-proximal) score. Position 1 = PAM-proximal.
    - If any mismatch in first 5 bases -> (0, True) meaning exclude.
    - Perfect in first 10 -> 50 pts.
    - First mismatch in positions 1..10 -> 50 - 5*position (min 0).
    - Mismatches only beyond 10 -> 50 pts.
    Returns (score, should_exclude).
    """
    positions = list(mismatch_positions_from_pam) if mismatch_positions_from_pam else []
    # Hard filter: any mismatch in first 5 bases -> exclude
    for p in positions:
        if 1 <= p <= SEED_FILTER_FIRST_N_BASES:
            return 0.0, True
    # Score: perfect in first 10 = 50; first mismatch in 1..10 reduces score
    first_mismatch_in_seed = None
    for p in positions:
        if 1 <= p <= 10:
            first_mismatch_in_seed = p
            break
    if first_mismatch_in_seed is None:
        return float(SEED_PERFECT_10_POINTS), False
    pts = SEED_PERFECT_10_POINTS - SEED_FIRST_MISMATCH_PENALTY_PER_POSITION * first_mismatch_in_seed
    return max(0.0, float(pts)), False


def _is_bd(rec: dict) -> bool:
    """On-target (0 mm) with any other *different* gene's TSS within BD_MAX_TSS_DISTANCE_BP.
    Same-gene isoforms (same gene_id) do not trigger BD. We check all TSSs within 2 kb,
    so even if the 2nd–10th nearest are same-gene isoforms, we still flag BD if an
    11th-nearest (or any) TSS from a different gene is within 2 kb.
    """
    if rec.get("mismatches", 0) != 0:
        return False
    return bool(rec.get("other_gene_tss_within_2kb"))


def _intended_gene_per_guide(records: list[dict]) -> dict[str, str]:
    """
    For each guide_id, determine intended gene:
    1. Parse guide_id as "GENE_N" (e.g. "SAXO1_0" → "SAXO1")
    2. Fallback: nearest_gene of first 0-mismatch hit
    """
    intended: dict[str, str] = {}
    # First pass: try to parse guide_id pattern "GENE_N" or "GENE_N_M"
    for rec in records:
        gid = rec.get("guide_id", "")
        if gid in intended:
            continue
        # Try pattern: "GENE_N" or "GENE_N_M" etc - extract GENE (part before first underscore)
        if "_" in gid:
            gene_from_id = gid.split("_")[0].strip()
            if gene_from_id:
                intended[gid] = gene_from_id
    # Second pass: for any guide_id not parsed, use first 0-mm hit's nearest_gene
    for rec in records:
        gid = rec.get("guide_id", "")
        if gid in intended:
            continue
        if rec.get("mismatches", 1) == 0:
            intended[gid] = (rec.get("nearest_gene") or "").strip()
    return intended


def score_and_filter_off_targets(
    records: list[dict],
    exclude_same_locus: bool = True,
) -> list[dict]:
    """
    Add score, flag (BD), same_locus; exclude records that fail the seed filter
    (any mismatch in first 5 bases from PAM).

    If exclude_same_locus is True (default), drop rows where the hit is in the
    same gene as the intended target (on-target / same-locus). That way the
    output CSV contains only off-targets.
    """
    intended = _intended_gene_per_guide(records)
    out = []
    for rec in records:
        rec = dict(rec)
        # Hard filter: mismatch in first 5 bases -> exclude
        positions = rec.get("mismatch_positions_from_pam") or []
        seed_pts, exclude = _seed_score(positions)
        if exclude:
            if os.environ.get("TARGETFINDER_DEBUG") and rec.get("guide_id") == "SAXO1_0" and str(rec.get("chr")) in ("9", "chr9"):
                print(f"[TargetFinder score] SAXO1_0 chr9 DROPPED by seed filter (mismatch in first 5 bases) positions={positions}", file=sys.stderr)
            continue

        tss_d = rec.get("tss_distance")
        tss_pts = _tss_score(tss_d)
        intended_gene = (intended.get(rec.get("guide_id", ""), "") or "").strip()
        nearest_gene = (rec.get("nearest_gene") or "").strip()
        # Compare gene names (case-insensitive for robustness)
        same_locus = bool(intended_gene and nearest_gene and intended_gene.upper() == nearest_gene.upper())

        if same_locus:
            rec["same_locus"] = True
            total = 0.0
            # BD = on-target with another gene's TSS within 2 kb (e.g. RRAGA near SAXO1). Keep those.
            bd = _is_bd(rec)
            if os.environ.get("TARGETFINDER_DEBUG") and rec.get("guide_id") == "SAXO1_0" and str(rec.get("chr")) in ("9", "chr9"):
                print(
                    f"[TargetFinder score] SAXO1_0 chr9 same_locus=True bd={bd} "
                    f"other_gene_tss_within_2kb={rec.get('other_gene_tss_within_2kb')} "
                    f"exclude_same_locus={exclude_same_locus} -> {'KEEP' if bd or not exclude_same_locus else 'DROP'}",
                    file=sys.stderr,
                )
            if exclude_same_locus and not bd:
                continue  # Drop same-locus only when it's not BD
        else:
            total = tss_pts + seed_pts
            rec["same_locus"] = False
            nm = rec.get("mismatches", 0)
            if isinstance(nm, int):
                total -= TOTAL_MISMATCH_PENALTY * nm
            total = max(0.0, total)
            bd = _is_bd(rec)

        if bd:
            rec["flag"] = "BD"
            rec["score"] = total + BD_BONUS
        else:
            rec["flag"] = ""
            rec["score"] = total

        out.append(rec)
    return out
