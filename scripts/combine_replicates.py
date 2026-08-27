#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
combine_replicates.py — Meta-combine v7.3 peak CSVs across biological replicates
==================================================================================

Takes K per-replicate peak CSVs (filterPeaks_normtoCTRL_v7.3.py outputs) and
emits a single combined CSV in which each unique cross-replicate site is one
row.  Combination uses:

  - **Stouffer's method** (sqrt-N weighted by default) to combine per-replicate
    p-values into a single meta-p, then BH-adjusted globally to a meta-q.
  - **Median and mean treated_sf** across the contributing replicates as
    effect-size summaries (median is the recommended primary effect-size
    metric; mean is included for reference).
  - **Replicate concordance counts**: how many replicates contributed any
    candidate at this site (`n_reps_in_meta`) and how many had `peak_pass=True`
    individually (`n_reps_passing`).

Cross-replicate sites are identified by the tuple
(`accession`, `tile_start`, `tile_end`, overlapping peak windows).  Peaks
from the same tile across replicates whose `peak_window_start` /
`peak_window_end` intervals overlap are clustered into a single site
(absorbs the ±1 nt biological fuzziness for the default w=3 window and
scales correctly if you ran v7.3 with a different `--pvalue-window`).
Peaks at the same genomic position but from **different tiles** are kept
as distinct sites: with v7.3's per-tile NB null, those peaks come from
independent fragment populations (different 3'-end pools) and have
test statistics calibrated to different tiles' parameters, so merging
them would conflate independent measurements.

The combined output carries `tile_start` and `tile_end` through to each
row so downstream tools (post_process_peaks.py, IGV) can identify which
tile each combined call came from.

Missing-replicate handling
--------------------------
A replicate's "missing" peak at a cross-replicate site means the position
never crossed v7.3's `min_adj_sf` seed threshold in that replicate, so no
candidate was emitted to that replicate's CSV.  How we treat this in the
meta-combination is a design choice (see `--missing-pvalue`):

  - `drop`     : exclude missing reps from Stouffer's; combine only
                 contributing reps' p-values. Statistically MNAR-biased
                 (positions are missing precisely *because* their statistic
                 was low, so dropping them inflates apparent significance).
                 Provided for backward compatibility.
  - `neutral`  : inject p = 0.5 (z = 0) for missing reps, with weight equal
                 to the median weight of the contributing reps. The missing
                 reps add to the weight denominator but contribute no signal,
                 so the combined p drifts toward less significant. **Default.**
                 The closest defensible compromise without going back to
                 the raw fragment data.
  - `one`      : inject p = 1 - 1e-15 for missing reps, with median weight.
                 The combined p plummets toward 1 — equivalent to poolr's
                 behaviour with a literal p=1 input. Use this only when you
                 want to require that every replicate explicitly produce a
                 candidate at the site.

Note that when v7.3 calls a peak with `peak_pass=False` (failed G1-G7
gates), the p-value was still computed and is a legitimate input to the
meta-combination. Such peaks DO contribute to `n_reps_found` and to
Stouffer's, but they do NOT contribute to `n_reps_passing`. The per-rep
fail reasons are surfaced in the `fail_reasons` summary column and in the
per-rep `peak_fail_reason_*` audit columns.

Use `--require-passing-in M` to require ≥ M replicates with `peak_pass=True`,
or `--require-found-in M` to require ≥ M replicates that produced any
candidate (passing or not).

Ranking
-------
The recommended sort is by ascending meta_qvalue with median_treated_sf as
a tie-breaker:

    sort_values(["meta_qvalue", "median_treated_sf"], ascending=[True, False])

This handles the floating-point underflow of very-significant peaks (which
all collapse to qvalue ≈ 0) by falling back to the effect-size column.

Usage
-----
    python scripts/combine_replicates.py R1.peaks.csv R2.peaks.csv R3.peaks.csv \\
        --out combined.csv \\
        --labels R1 R2 R3

Pass `--fasta ref.fa` to also append downstream_seq* columns (see --fasta
help) to the combined CSV directly, reusing extract_downstream_seqs.py's
window-extraction logic.
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats

from extract_downstream_seqs import extract_window, read_fasta


# ─────────────────────────────────────────────────────────────────────────────
# I/O
# ─────────────────────────────────────────────────────────────────────────────

REQUIRED_COLS = {
    "accession", "tile_id", "tile_start", "tile_end",
    "stall_position", "stall_site_genomic_position",
    "peak_window_start", "peak_window_end",
    "treated_sf", "tile_reads_treated", "pvalue", "qvalue", "peak_pass",
}


def read_peak_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    missing = REQUIRED_COLS - set(df.columns)
    if missing:
        raise SystemExit(
            f"[combine] {path} is missing columns: {sorted(missing)}.  "
            "Is this a filterPeaks_normtoCTRL_v7.3.py output?"
        )
    df["accession"] = df["accession"].astype(str)
    for col in ("tile_id", "tile_start", "tile_end", "stall_position",
                "stall_site_genomic_position",
                "peak_window_start", "peak_window_end", "tile_reads_treated"):
        df[col] = df[col].astype(np.int64)
    df["peak_pass"] = df["peak_pass"].astype(bool)
    return df


# ─────────────────────────────────────────────────────────────────────────────
# Position alignment by window overlap
# ─────────────────────────────────────────────────────────────────────────────


def cluster_overlapping(peaks_df: pd.DataFrame) -> List[List[int]]:
    """
    Greedy clustering of peaks within one accession by overlapping
    [peak_window_start, peak_window_end] intervals.

    `peaks_df` must already be sorted by peak_window_start.
    Returns a list of clusters, each a list of integer row indices into
    peaks_df.
    """
    clusters: List[List[int]] = []
    if len(peaks_df) == 0:
        return clusters
    starts = peaks_df["peak_window_start"].to_numpy(np.int64)
    ends = peaks_df["peak_window_end"].to_numpy(np.int64)

    current = [0]
    current_end = int(ends[0])
    for i in range(1, len(peaks_df)):
        if int(starts[i]) <= current_end:
            current.append(i)
            current_end = max(current_end, int(ends[i]))
        else:
            clusters.append(current)
            current = [i]
            current_end = int(ends[i])
    clusters.append(current)
    return clusters


# ─────────────────────────────────────────────────────────────────────────────
# Statistical combination
# ─────────────────────────────────────────────────────────────────────────────


def stouffer_combine(pvals: np.ndarray, weights: np.ndarray) -> Tuple[float, float]:
    """
    Stouffer's weighted combination of p-values.

    Returns (combined_pvalue, combined_z).  Equivalent to combining via
    z = norm.isf(p) (numerically stable in the deep tail unlike ppf(1-p)),
    weighted-summed, then norm.sf back to a p-value.

    p-values are clipped to [1e-300, 1-1e-15] to keep z finite.  With one
    replicate the result is the input p (Stouffer's collapses correctly).
    """
    pvals = np.clip(np.asarray(pvals, dtype=np.float64), 1e-300, 1.0 - 1e-15)
    weights = np.asarray(weights, dtype=np.float64)
    if pvals.size == 0:
        return float("nan"), float("nan")
    z = stats.norm.isf(pvals)
    z_combined = float(np.sum(weights * z) / np.sqrt(np.sum(weights ** 2)))
    return float(stats.norm.sf(z_combined)), z_combined


def bh_qvalues(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values, monotone-non-increasing from the right."""
    pvals = np.asarray(pvals, dtype=np.float64)
    n = pvals.size
    if n == 0:
        return pvals
    order = np.argsort(pvals)
    ranked = pvals[order]
    raw = ranked * n / (np.arange(n) + 1)
    adj_sorted = np.minimum.accumulate(raw[::-1])[::-1]
    adj_sorted = np.clip(adj_sorted, 0.0, 1.0)
    out = np.empty(n, dtype=np.float64)
    out[order] = adj_sorted
    return out


# ─────────────────────────────────────────────────────────────────────────────
# Per-cluster aggregation
# ─────────────────────────────────────────────────────────────────────────────


# Per-replicate columns to surface in the output (in addition to whichever
# columns drive the meta-statistics).  Kept narrow on purpose to keep the
# combined table readable; users can join back to the per-rep CSVs for more.
PER_REP_COLS = (
    "pvalue",
    "qvalue",
    "treated_sf",
    "adj_sf",
    "tile_reads_treated",
    "stall_reads_treated",
    "peak_pass",
    "peak_fail_reason",
)


def aggregate_cluster(
    accession: str,
    tile_start: int,
    tile_end: int,
    cluster_rows: List[Tuple[str, pd.Series]],
    all_labels: List[str],
    weight_scheme: str,
    missing_pvalue: str,
) -> Tuple[dict, Dict[str, pd.Series]]:
    """
    Build one output row from a cluster of (label, peak_row) tuples.

    For each contributing replicate, take the most-significant peak in the
    cluster (smallest pvalue) — guards against the rare case that the
    same replicate has multiple within-tolerance candidates.

    `all_labels` is the full list of replicate labels (in the user-supplied
    order) so we know which replicates are MISSING from this cluster.  The
    `missing_pvalue` policy decides what to inject for those:

      - "drop"    : exclude missing reps from Stouffer's (n_used = n_found)
      - "neutral" : inject p = 0.5 (z = 0) and a median-weight; missing reps
                    contribute no signal but inflate the weight denominator,
                    so the combined p drifts toward less significant
      - "one"     : inject p = 1 - 1e-15 (z very negative) — combined p
                    plummets toward 1 (matches poolr behaviour with p=1)
    """
    # Per-replicate selection: best (smallest p) peak per label.
    best_per_rep: Dict[str, pd.Series] = {}
    for label, row in cluster_rows:
        if label not in best_per_rep or row["pvalue"] < best_per_rep[label]["pvalue"]:
            best_per_rep[label] = row

    # Cluster-level stats from the contributing replicates only.
    contribs = list(best_per_rep.items())
    contrib_pvals = np.array([row["pvalue"] for _, row in contribs], dtype=np.float64)
    contrib_Ns = np.array([row["tile_reads_treated"] for _, row in contribs], dtype=np.float64)
    treated_sfs = np.array([row["treated_sf"] for _, row in contribs], dtype=np.float64)
    positions = np.array([row["stall_position"] for _, row in contribs], dtype=np.int64)
    window_starts = np.array([row["peak_window_start"] for _, row in contribs], dtype=np.int64)
    window_ends = np.array([row["peak_window_end"] for _, row in contribs], dtype=np.int64)
    n_passing = int(sum(bool(row["peak_pass"]) for _, row in contribs))
    n_found = len(contribs)

    if weight_scheme == "sqrt-n":
        contrib_weights = np.sqrt(contrib_Ns)
    elif weight_scheme == "uniform":
        contrib_weights = np.ones_like(contrib_Ns)
    else:
        raise SystemExit(f"[combine] Unknown --weight: {weight_scheme!r}")

    # ── Apply missing-pvalue policy ─────────────────────────────────────────
    missing_labels = [lbl for lbl in all_labels if lbl not in best_per_rep]
    if missing_pvalue == "drop" or len(missing_labels) == 0:
        meta_pvals = contrib_pvals
        meta_weights = contrib_weights
    else:
        if missing_pvalue == "neutral":
            inject_p = 0.5
        elif missing_pvalue == "one":
            inject_p = 1.0 - 1e-15
        else:
            raise SystemExit(f"[combine] Unknown --missing-pvalue: {missing_pvalue!r}")
        # Use the median of contributing weights as the weight for missing
        # reps — keeps the missing reps from dominating either way and is
        # the most defensible heuristic when their N is unknown.
        inject_w = float(np.median(contrib_weights)) if contrib_weights.size > 0 else 1.0
        n_missing = len(missing_labels)
        meta_pvals = np.concatenate([contrib_pvals, np.full(n_missing, inject_p)])
        meta_weights = np.concatenate([contrib_weights, np.full(n_missing, inject_w)])

    n_used = int(meta_pvals.size)
    meta_p, meta_z = stouffer_combine(meta_pvals, meta_weights)

    # Representative position: median across contributing peaks, rounded to int.
    rep_pos = int(np.round(np.median(positions)))
    # All contributing rows share the same accession (and thus the same parent-
    # genome offset), so derive offset from any one of them and apply to rep_pos.
    _ref = contribs[0][1]
    genome_offset = int(_ref["stall_site_genomic_position"]) - int(_ref["stall_position"])
    rep_pos_genomic = rep_pos + genome_offset
    # Cluster span: union of all per-rep windows.
    cluster_window_start = int(window_starts.min())
    cluster_window_end = int(window_ends.max())

    # ── fail_reasons summary (per-rep label : fail string) ─────────────────
    fail_parts: List[str] = []
    for lbl in all_labels:
        if lbl in best_per_rep:
            row = best_per_rep[lbl]
            if not bool(row["peak_pass"]):
                reason = row.get("peak_fail_reason", "")
                reason = "" if (reason is None or (isinstance(reason, float) and np.isnan(reason))) else str(reason)
                fail_parts.append(f"{lbl}:{reason if reason else 'failed'}")
        else:
            fail_parts.append(f"{lbl}:not_called")
    fail_reasons = "; ".join(fail_parts)

    out = {
        # Identifier columns leading the row so they're easy to spot when
        # the CSV is opened in a spreadsheet / used as a join key downstream.
        #   peak_position_id  — unique per (accession, stall_position); the
        #                       canonical name for a cross-rep stall site.
        #   peak_id           — unique per row (adds tile span); useful when
        #                       multiple overlapping tiles call the same
        #                       position and you need to distinguish which
        #                       tile's evidence produced this row.
        "peak_id": f"{accession}|{rep_pos}|{tile_start}-{tile_end}",
        "peak_position_id": f"{accession}|{rep_pos}",
        "peak_genomic_position_id": f"{accession}|{rep_pos_genomic}",
        "accession": accession,
        "tile_id": int(contribs[0][1]["tile_id"]),
        "tile_start": tile_start,
        "tile_end": tile_end,
        "stall_position": rep_pos,
        "stall_site_genomic_position": rep_pos_genomic,
        # stall_site_rank is recomputed across the combined candidate set
        # after build_combined; the per-rep value is not meaningful here.
        "cluster_window_start": cluster_window_start,
        "cluster_window_end": cluster_window_end,
        "meta_pvalue": meta_p,
        "meta_z": meta_z,
        "n_reps_found": n_found,
        "n_reps_passing": n_passing,
        "n_reps_used_in_meta": n_used,
        "median_treated_sf": float(np.median(treated_sfs)),
        "mean_treated_sf": float(np.mean(treated_sfs)),
        "median_tile_reads_treated": float(np.median(contrib_Ns)),
        "fail_reasons": fail_reasons,
        # meta_qvalue gets filled in after we have the full table (BH is global).
    }
    return out, best_per_rep


def build_combined(
    per_rep: Dict[str, pd.DataFrame],
    weight_scheme: str,
    missing_pvalue: str,
    min_reps_passing_for_rank: int,
    max_ranked_per_accession: int,
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Combine per-replicate peak DataFrames into one row per cross-rep site.

    Returns (combined_df, ordered_labels).  Per-rep audit columns are added
    in label order so the output is deterministic.
    """
    labels = list(per_rep.keys())

    # Tag each row with its replicate label and pool everything per
    # (accession, tile_start, tile_end) — the tile is part of the site
    # identity because peaks in different tiles at the same genomic position
    # come from different fragment populations (different 3'-end pools)
    # and have per-tile-calibrated test statistics.  Merging them across
    # tiles would double-count and conflate independent measurements.
    per_tile: Dict[Tuple[str, int, int], List[Tuple[str, pd.Series]]] = defaultdict(list)
    for label, df in per_rep.items():
        for _, row in df.iterrows():
            key = (str(row["accession"]), int(row["tile_start"]), int(row["tile_end"]))
            per_tile[key].append((label, row))

    combined_rows: List[dict] = []
    per_rep_lookups: List[Dict[str, pd.Series]] = []

    for (acc, tile_start, tile_end) in sorted(per_tile.keys()):
        # Within one tile, cluster peaks across replicates by overlapping
        # peak_window_start/peak_window_end intervals (handles ±1 nt fuzziness).
        labels_local: List[str] = []
        rows_local: List[pd.Series] = []
        for label, row in per_tile[(acc, tile_start, tile_end)]:
            labels_local.append(label)
            rows_local.append(row)
        small = pd.DataFrame({
            "_label": labels_local,
            "peak_window_start": [r["peak_window_start"] for r in rows_local],
            "peak_window_end":   [r["peak_window_end"]   for r in rows_local],
            "_idx": list(range(len(rows_local))),
        }).sort_values("peak_window_start", kind="mergesort").reset_index(drop=True)

        clusters = cluster_overlapping(small)

        for cluster in clusters:
            cluster_rows = [(small["_label"].iloc[i], rows_local[int(small["_idx"].iloc[i])])
                            for i in cluster]
            row, best_per_rep = aggregate_cluster(
                acc, tile_start, tile_end, cluster_rows,
                all_labels=labels,
                weight_scheme=weight_scheme,
                missing_pvalue=missing_pvalue,
            )
            combined_rows.append(row)
            per_rep_lookups.append(best_per_rep)

    df = pd.DataFrame(combined_rows)

    # Global BH on meta_pvalue → meta_qvalue.  Insert immediately after
    # meta_pvalue so the two are visually adjacent in the output (rather
    # than meta_qvalue being appended at the end).
    if not df.empty:
        meta_q = bh_qvalues(df["meta_pvalue"].to_numpy())
        meta_p_idx = df.columns.get_loc("meta_pvalue")
        df.insert(meta_p_idx + 1, "meta_qvalue", meta_q)

        # Recompute stall_site_rank with eligibility + cap:
        #   - Eligible peaks: n_reps_passing >= min_reps_passing_for_rank
        #     (i.e. "not failed" in at least M replicates).
        #   - Cap per accession: at most max_ranked_per_accession peaks get a
        #     numeric rank, chosen as the top-N by meta_pvalue (most
        #     significant). Ties broken by smaller stall_site_genomic_position
        #     (more 5') so the choice is deterministic.
        #   - Within the chosen subset, rank by stall_site_genomic_position
        #     ascending (5'-most = 1), dense — tied positions share a rank.
        #   - Eligible but not in the top-N → "-".
        #   - Not eligible (n_reps_passing < M) → NaN.
        # Column dtype is object/string because of the "-" sentinel.
        ranks = pd.Series([pd.NA] * len(df), index=df.index, dtype="object")
        eligible_mask = df["n_reps_passing"] >= min_reps_passing_for_rank
        for acc, g_idx in df[eligible_mask].groupby("accession", sort=False).groups.items():
            g = df.loc[g_idx]
            chosen = (
                g.sort_values(
                    ["meta_pvalue", "stall_site_genomic_position"],
                    ascending=[True, True],
                    kind="mergesort",
                )
                .head(max_ranked_per_accession)
            )
            ranks.loc[g.index.difference(chosen.index)] = "-"
            chosen_ranks = (
                chosen["stall_site_genomic_position"]
                .rank(method="dense", ascending=True)
                .astype(int)
            )
            for idx, r in chosen_ranks.items():
                ranks.loc[idx] = str(int(r))
        df["stall_site_rank"] = ranks

    # Append per-rep audit columns.
    for col in PER_REP_COLS:
        for label in labels:
            out_col = f"{col}_{label}"
            df[out_col] = [
                (lookup[label][col] if label in lookup else np.nan)
                for lookup in per_rep_lookups
            ]

    # Surface the per-rep tile_reads_treated_* values immediately after the
    # SF / read-count summary block — coverage variability across reps is
    # the most decision-relevant per-rep audit info ("did one rep have way
    # fewer reads at this tile?") and it should be visible without scrolling.
    # Anchor on median_tile_reads_treated if present (newer output), else
    # mean_treated_sf, so the summary stats stay grouped together before the
    # per-rep details.
    if not df.empty:
        reads_cols = [f"tile_reads_treated_{lbl}" for lbl in labels
                      if f"tile_reads_treated_{lbl}" in df.columns]
        if reads_cols:
            for anchor in ("median_tile_reads_treated", "mean_treated_sf"):
                if anchor in df.columns:
                    anchor_idx = df.columns.get_loc(anchor) + 1
                    other_cols = [c for c in df.columns if c not in reads_cols]
                    new_order = other_cols[:anchor_idx] + reads_cols + other_cols[anchor_idx:]
                    df = df[new_order]
                    break

    return df, labels


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "peaks_csvs", nargs="+", type=Path,
        help="K per-replicate v7.3 peak CSVs (positional).",
    )
    parser.add_argument(
        "--out", required=True, type=Path,
        help="Output combined CSV path.",
    )
    parser.add_argument(
        "--labels", nargs="+", default=None,
        help="One label per replicate (used as the suffix on per-rep audit columns "
        "and in summary output). Defaults to file stems.",
    )
    parser.add_argument(
        "--weight", choices=["sqrt-n", "uniform"], default="sqrt-n",
        help="Stouffer's weighting scheme: sqrt-n (Liptak; weights by sqrt(tile_reads_treated)) "
        "or uniform. Default sqrt-n.",
    )
    parser.add_argument(
        "--missing-pvalue", choices=["drop", "neutral", "one"], default="neutral",
        help="Policy for replicates that did not produce a candidate at a cross-rep "
        "site (i.e. the position never crossed v7.3's --min-adj-sf seed threshold "
        "in that replicate). 'drop' (legacy) excludes them from Stouffer's; "
        "'neutral' injects p=0.5 (z=0) so they dilute the combined signal without "
        "killing it (default; the right compromise without going back to fragments); "
        "'one' injects p≈1, which forces the combined p toward 1 — equivalent to "
        "poolr::stouffer with a literal p=1 input. See the script docstring for the "
        "rationale; ideally you want 'fragments' mode (per-position raw recompute) "
        "but that's not yet implemented.",
    )
    parser.add_argument(
        "--require-passing-in", type=int, default=0,
        help="Drop combined sites where fewer than M replicates had peak_pass=True. "
        "Default 0 (no filter; n_reps_passing is annotated regardless).",
    )
    parser.add_argument(
        "--require-found-in", type=int, default=0,
        help="Drop combined sites where fewer than M replicates produced any candidate "
        "at the position. Default 0 (no filter; n_reps_found is annotated regardless).",
    )
    parser.add_argument(
        "--max-meta-qvalue", type=float, default=None,
        help="Optional: drop combined sites with meta_qvalue above this threshold.",
    )
    parser.add_argument(
        "--min-reps-passing-for-rank", type=int, default=2,
        help="Minimum n_reps_passing for a peak to be eligible for a numeric "
        "stall_site_rank. Peaks below this stay in the output but get stall_site_rank=NaN. "
        "Default 2 (e.g. 2-of-3 reps must pass).",
    )
    parser.add_argument(
        "--max-ranked-per-accession", type=int, default=5,
        help="Maximum number of peaks per accession that receive a numeric "
        "stall_site_rank. When more eligible peaks exist, the top-N by meta_pvalue "
        "(most significant) are kept; the rest are marked '-'. Within the chosen "
        "subset the rank itself is by stall_site_genomic_position (5'-most = 1). "
        "Default 5.",
    )
    parser.add_argument(
        "--sort", choices=["meta_qvalue", "meta_pvalue", "median_treated_sf", "none"],
        default="meta_qvalue",
        help="Primary sort key for the output. With meta_qvalue (default), median_treated_sf "
        "is used as a tie-breaker (descending) to handle floating-point underflow at the "
        "very-significant end.",
    )
    parser.add_argument(
        "--fasta", type=Path, default=None,
        help="Optional reference FASTA. When given, the downstream sequence "
        "around each combined site's stall_position is extracted (reusing "
        "extract_downstream_seqs.py's window logic) and appended as "
        "downstream_seq / downstream_seq_start / downstream_seq_end / "
        "downstream_stall_at_offset / downstream_seq_note columns — one "
        "sequence per output row, unfiltered (combine_replicates.py already "
        "emits one row per cross-rep site, so there's no separate tiebreak "
        "or dedup step here). downstream_seq_note is 'ok', "
        "'missing_accession', or 'empty_window'.",
    )
    parser.add_argument(
        "--upstream", type=int, default=5,
        help="Bases of upstream context before stall_position when --fasta is given. "
        "Clipped to tile_start. Default 5.",
    )
    parser.add_argument(
        "--max-length", type=int, default=80,
        help="Max downstream length in nt from stall_position when --fasta is given. "
        "Clipped to tile_end. Default 80.",
    )
    args = parser.parse_args()

    if args.labels is not None and len(args.labels) != len(args.peaks_csvs):
        raise SystemExit(
            f"[combine] --labels has {len(args.labels)} entries but "
            f"{len(args.peaks_csvs)} peak CSVs were supplied; counts must match."
        )

    labels = args.labels or [p.stem for p in args.peaks_csvs]
    if len(set(labels)) != len(labels):
        raise SystemExit(f"[combine] Replicate labels must be unique; got {labels}")

    print(f"[combine] Reading {len(args.peaks_csvs)} replicate(s)...", file=sys.stderr)
    per_rep: Dict[str, pd.DataFrame] = {}
    for label, path in zip(labels, args.peaks_csvs):
        df = read_peak_csv(path)
        per_rep[label] = df
        print(f"[combine]   {label}: {len(df)} candidate peaks ({int(df['peak_pass'].sum())} passing)  <- {path}",
              file=sys.stderr)

    print(
        f"[combine] Clustering across replicates and combining "
        f"(weight={args.weight}, missing-pvalue={args.missing_pvalue})...",
        file=sys.stderr,
    )
    combined, ordered_labels = build_combined(
        per_rep,
        weight_scheme=args.weight,
        missing_pvalue=args.missing_pvalue,
        min_reps_passing_for_rank=args.min_reps_passing_for_rank,
        max_ranked_per_accession=args.max_ranked_per_accession,
    )

    if combined.empty:
        print("[combine] No cross-replicate sites produced.", file=sys.stderr)
        combined.to_csv(args.out, index=False)
        return

    # Summary distribution of replicate participation, before any filtering.
    found_counts = combined["n_reps_found"].value_counts().sort_index()
    pass_counts = combined["n_reps_passing"].value_counts().sort_index()
    print(
        "[combine] Cross-rep site counts by n_reps_found: "
        + ", ".join(f"{int(k)}={int(v)}" for k, v in found_counts.items()),
        file=sys.stderr,
    )
    print(
        "[combine] Cross-rep site counts by n_reps_passing: "
        + ", ".join(f"{int(k)}={int(v)}" for k, v in pass_counts.items()),
        file=sys.stderr,
    )

    # Apply optional filters.
    n_before = len(combined)
    if args.require_passing_in > 0:
        combined = combined[combined["n_reps_passing"] >= args.require_passing_in]
    if args.require_found_in > 0:
        combined = combined[combined["n_reps_found"] >= args.require_found_in]
    if args.max_meta_qvalue is not None:
        combined = combined[combined["meta_qvalue"] <= args.max_meta_qvalue]
    n_after = len(combined)
    if n_after < n_before:
        print(
            f"[combine] Filters reduced {n_before} -> {n_after} sites "
            f"(--require-passing-in={args.require_passing_in}, "
            f"--require-found-in={args.require_found_in}, "
            f"--max-meta-qvalue={args.max_meta_qvalue}).",
            file=sys.stderr,
        )

    # Sort.
    if args.sort == "meta_qvalue":
        combined = combined.sort_values(
            ["meta_qvalue", "median_treated_sf"], ascending=[True, False], kind="mergesort"
        )
    elif args.sort == "meta_pvalue":
        combined = combined.sort_values(
            ["meta_pvalue", "median_treated_sf"], ascending=[True, False], kind="mergesort"
        )
    elif args.sort == "median_treated_sf":
        combined = combined.sort_values("median_treated_sf", ascending=False, kind="mergesort")
    # "none": leave in build order (accession-then-window-position).

    if args.fasta is not None:
        print(f"[combine] Reading FASTA: {args.fasta}", file=sys.stderr)
        fasta = read_fasta(args.fasta)

        def extract_row(row: pd.Series) -> pd.Series:
            acc = str(row["accession"])
            if acc not in fasta:
                return pd.Series({"downstream_seq": "", "downstream_seq_start": pd.NA,
                                   "downstream_seq_end": pd.NA, "downstream_stall_at_offset": pd.NA,
                                   "downstream_seq_note": "missing_accession"})
            stall = int(row["stall_position"])
            t_start = int(row["tile_start"])
            t_end = int(row["tile_end"])
            seq_start, seq_end, sequence = extract_window(
                fasta[acc], stall, t_start, t_end,
                upstream=args.upstream, max_length=args.max_length,
            )
            if not sequence:
                return pd.Series({"downstream_seq": "", "downstream_seq_start": seq_start,
                                   "downstream_seq_end": seq_end, "downstream_stall_at_offset": pd.NA,
                                   "downstream_seq_note": "empty_window"})
            return pd.Series({"downstream_seq": sequence, "downstream_seq_start": seq_start,
                               "downstream_seq_end": seq_end,
                               "downstream_stall_at_offset": stall - seq_start,
                               "downstream_seq_note": "ok"})

        seq_cols = combined.apply(extract_row, axis=1)
        combined = pd.concat([combined, seq_cols], axis=1)
        note_counts = combined["downstream_seq_note"].value_counts().to_dict()
        print(f"[combine] Downstream sequence extraction: {note_counts}", file=sys.stderr)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(args.out, index=False)
    print(
        f"[combine] Wrote {len(combined)} cross-rep site(s)  ->  {args.out}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
