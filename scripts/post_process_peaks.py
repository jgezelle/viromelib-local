#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
post_process_peaks.py — Filter, annotate, and (optionally) collapse peaks
================================================================================

Accepts either input format and dispatches automatically:

  - **Per-replicate** v7.3 peak CSV (filterPeaks_normtoCTRL_v7.3.py output) —
    filters dispatch to `peak_pass`, `treated_sf`, `tile_reads_treated`.
  - **Cross-replicate combined** CSV (combine_replicates.py output, identified
    by the presence of the `meta_pvalue` column) — filters dispatch to
    `n_reps_passing`, `median_treated_sf`, and the median of the per-replicate
    `tile_reads_treated_*` columns.

Performs three post-processing steps:

  1. FILTER. Keep only rows with peak_pass == True, then apply the optional
     --min-treated-sf (on treated_sf) and --min-reads (on tile_reads_treated)
     thresholds.  Both default to 0, i.e. peak_pass alone is the only filter
     unless overridden.  These filters apply to both outputs.

  2. ANNOTATE. Add taxonomy columns (`family`, `genus`, `species` — whichever
     are present) by looking up each peak's accession in a taxonomy CSV
     (e.g. virus_stats.csv).  Peak accessions like
     "NC_076379.1|2658-2844_inferred_3UTR" are split on "|" and matched on the
     leading canonical NCBI accession.  `family` is required in the lookup CSV;
     `genus` and `species` are added automatically when present.  Optionally
     also annotate Infernal hits from one or more .tbl files (4 diagnostic
     columns per tbl + any_hit; same scheme as the legacy
     annotate_infernal_hits.py).

  3. EMIT. Write a per-peak "raw" CSV with all annotations attached, and
     optionally a per-accession "collapsed" CSV in which every per-peak field
     is semicolon-joined within the row.

Inputs
------
  PEAKS_CSV               positional — output of filterPeaks_normtoCTRL_v7.3.py
  --family-csv PATH       required — must have at least 'accession' and 'family' columns
  --out-raw PATH          required — per-peak output
  --out-collapsed PATH    optional — per-accession output
  --infernal-tbl PATH     repeatable — Infernal .tbl files (uses peak accession verbatim)
  --min-treated-sf FLOAT  optional (default 0) — drop peaks with treated_sf below this
  --min-reads INT         optional (default 0) — drop peaks with tile_reads_treated below this
  --additional-columns C1,C2,...
                          extra peak columns to also semicolon-join into the
                          collapsed output

Collapsed output schema
-----------------------
One row per peak `accession`. Single-valued fields (one per accession): the
full peak accession itself, family, peak_count.  Per-peak fields (semicolon-
joined across the accession's peaks, in row order): stall_position, adj_sf,
tile_reads_treated, pvalue, qvalue, any_hit (if --infernal-tbl was used),
plus anything supplied via --additional-columns.
"""

import functools
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import typer

app = typer.Typer(add_completion=False)

INPUT_FILE_OPTS = {
    "exists": True,
    "file_okay": True,
    "dir_okay": False,
    "readable": True,
    "resolve_path": True,
}
OUTPUT_FILE_OPTS = {
    "file_okay": True,
    "dir_okay": False,
    "writable": True,
    "resolve_path": True,
}
InputFileArg = functools.partial(typer.Argument, **INPUT_FILE_OPTS)
InputFileOpt = functools.partial(typer.Option, **INPUT_FILE_OPTS)
OutputFileOpt = functools.partial(typer.Option, **OUTPUT_FILE_OPTS, default=None)


# ─────────────────────────────────────────────────────────────────────────────
# Family lookup (virus_stats.csv style)
# ─────────────────────────────────────────────────────────────────────────────


# Taxonomy columns to pull from --family-csv if present.  `family` is required;
# `genus` and `species` are added when the lookup CSV provides them, otherwise
# silently skipped.  Order here is the order they appear in the output (left
# to right).
TAXONOMY_COLUMNS = ("family", "genus", "species")


def add_family(df: pd.DataFrame, family_csv: Path) -> pd.DataFrame:
    """
    Add taxonomy columns (`family`, `genus`, `species` — whichever are present
    in the lookup CSV) to df by looking up the canonical NCBI accession
    (peak accession split on '|', first part) in `family_csv`.

    The lookup CSV must have at least 'accession' and 'family' columns;
    'genus' and 'species' are added if present and silently skipped if not.
    Non-identical duplicates per accession trigger a warning and the first
    row wins.  Peaks whose canonical accession is not in the lookup get
    NaN for every taxonomy column, with a single summary warning.
    """
    fam = pd.read_csv(family_csv, dtype=str)
    if "accession" not in fam.columns:
        raise SystemExit(f"[postproc] --family-csv {family_csv} missing 'accession' column.")
    if "family" not in fam.columns:
        raise SystemExit(f"[postproc] --family-csv {family_csv} missing 'family' column.")

    available_tax = [c for c in TAXONOMY_COLUMNS if c in fam.columns]
    skipped_tax = [c for c in TAXONOMY_COLUMNS if c not in fam.columns]
    if skipped_tax:
        print(
            f"[postproc] Note: {family_csv.name} does not contain {skipped_tax}; only {available_tax} will be added.",
            file=sys.stderr,
        )

    # Detect non-identical duplicates per accession across any taxonomy column.
    fam_first = fam.drop_duplicates("accession", keep="first")
    for col in available_tax:
        dup_per_acc = fam.groupby("accession")[col].nunique()
        conflicting = dup_per_acc[dup_per_acc > 1]
        if len(conflicting) > 0:
            print(
                f"[postproc] WARN: {len(conflicting)} accession(s) in {family_csv.name} have "
                f"non-identical '{col}' entries; first occurrence wins. "
                f"Examples: {conflicting.head(3).index.tolist()}",
                file=sys.stderr,
            )

    canonical = df["accession"].astype(str).str.split("|", n=1).str[0]

    df = df.copy()
    for col in available_tax:
        col_lookup = dict(zip(fam_first["accession"], fam_first[col]))
        df[col] = canonical.map(col_lookup).values

    # One summary warning, based on family (presence/absence of family is the
    # canonical "did the accession exist in the lookup at all" check).
    n_missing = int(df["family"].isna().sum())
    if n_missing > 0:
        unique_missing = canonical[df["family"].isna()].unique()
        examples = list(unique_missing[:3])
        print(
            f"[postproc] WARN: {n_missing} peak(s) across {len(unique_missing)} unique "
            f"accession(s) have no entry in {family_csv.name}; "
            f"{', '.join(available_tax)} = NaN. "
            f"Examples: {examples}",
            file=sys.stderr,
        )

    return df


# ─────────────────────────────────────────────────────────────────────────────
# Infernal annotation
# Behaviour copy of annotate_infernal_hits.py: 4 diagnostic columns per tbl
# (best E-value, hit start/end relative to apex, tile-truncation length) and
# a single any_hit boolean across all tbls.  Matching uses the peak accession
# verbatim (i.e. the same string the .tbl file uses as its target).
# ─────────────────────────────────────────────────────────────────────────────


def parse_tbl(tbl_path: Path, evalue_threshold: float) -> Dict[str, List[tuple]]:
    """Parse Infernal .tbl into accession -> list of (seq_from_0based, seq_to_0based, evalue)."""
    hits: Dict[str, List[tuple]] = defaultdict(list)
    with open(tbl_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 16:
                continue
            target = parts[0]
            strand = parts[9]
            seq_from = int(parts[7])
            seq_to = int(parts[8])
            try:
                evalue = float(parts[15])
            except ValueError:
                continue
            if strand != "+":
                continue
            if evalue > evalue_threshold:
                continue
            hits[target].append((seq_from - 1, seq_to - 1, evalue))
    return hits


def _tbl_stem(path: Path) -> str:
    return os.path.splitext(os.path.basename(str(path)))[0]


def annotate_infernal(
    df: pd.DataFrame,
    tbl_paths: List[Path],
    evalue: float,
    downstream: int,
    min_downstream_frac: float,
) -> pd.DataFrame:
    """For each tbl file, attach 4 columns + a single `any_hit` across all tbls."""
    df = df.copy()

    for tbl_path in tbl_paths:
        stem = _tbl_stem(tbl_path)
        hits_by_acc = parse_tbl(tbl_path, evalue)

        hit_col = f"{stem}_hit"
        start_col = f"{stem}_start"
        end_col = f"{stem}_end"
        trunc_col = f"{stem}_tile_truncated"

        best_evalue, hit_starts, hit_ends, trunc_length = [], [], [], []

        for _, row in df.iterrows():
            accession = row["accession"]
            stall = int(row["stall_position"])
            tile_start = int(row["tile_start"])
            tile_end = int(row["tile_end"])

            candidates = [
                (sf, st, ev)
                for sf, st, ev in hits_by_acc.get(accession, [])
                if sf <= stall + downstream
                and st >= stall
                and max(0, st - stall) / (st - sf + 1) >= min_downstream_frac
            ]

            if not candidates:
                best_evalue.append(float("nan"))
                hit_starts.append(float("nan"))
                hit_ends.append(float("nan"))
                trunc_length.append(0)
            else:
                sf, st, ev = min(candidates, key=lambda x: x[2])
                best_evalue.append(ev)
                hit_starts.append(sf - stall)
                hit_ends.append(st - stall)
                left_overrun = max(0, tile_start - sf)
                right_overrun = max(0, st - (tile_end - 1))
                trunc_length.append(left_overrun + right_overrun)

        df[hit_col] = best_evalue
        df[start_col] = hit_starts
        df[end_col] = hit_ends
        df[trunc_col] = trunc_length

    hit_cols = [f"{_tbl_stem(p)}_hit" for p in tbl_paths]
    df["any_hit"] = df[hit_cols].notna().any(axis=1) if hit_cols else False
    return df


# ─────────────────────────────────────────────────────────────────────────────
# Per-accession collapse
# ─────────────────────────────────────────────────────────────────────────────


# Columns reduced by OR across an accession's peaks rather than semicolon-joined.
# These are inherently boolean per-peak summaries where the question for the
# accession is "does ANY peak satisfy it?" — `any_hit` being the canonical case.
OR_REDUCE_COLUMNS = {"any_hit"}


def collapse_by_accession(
    df: pd.DataFrame,
    join_cols: List[str],
) -> pd.DataFrame:
    """
    One row per accession.  Single-valued columns (accession, family) appear
    once.  Columns in `join_cols` that are listed in OR_REDUCE_COLUMNS are
    reduced with `any()` (boolean OR across peaks); everything else is
    semicolon-joined in row order.  Adds `peak_count`.
    """
    if df.empty:
        leading = [c for c in TAXONOMY_COLUMNS if c in df.columns]
        cols = leading + ["accession", "peak_count"] + [c for c in join_cols if c in df.columns]
        return pd.DataFrame(columns=cols)

    def _join(series: pd.Series) -> str:
        return ";".join("" if pd.isna(v) else str(v) for v in series.tolist())

    def _or_reduce(series: pd.Series) -> bool:
        return bool(series.fillna(False).astype(bool).any())

    grouped = df.groupby("accession", sort=False)
    rows = []
    for acc, g in grouped:
        # Taxonomy columns are single-valued per accession (every peak in the
        # accession looked up the same row from --family-csv), so emit each
        # once if present in the input.  Family/genus/species ordering matches
        # TAXONOMY_COLUMNS.
        row = {}
        for tax_col in TAXONOMY_COLUMNS:
            if tax_col in g.columns:
                row[tax_col] = g[tax_col].iloc[0]
        row["accession"] = acc
        row["peak_count"] = len(g)
        for col in join_cols:
            if col not in g.columns:
                continue
            if col in OR_REDUCE_COLUMNS:
                row[col] = _or_reduce(g[col])
            else:
                row[col] = _join(g[col])
        rows.append(row)
    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────


@app.command()
def main(
    peaks_csv: Path = InputFileArg(help="Peak CSV from filterPeaks_normtoCTRL_v7.3.py."),
    out_raw: Path = typer.Option(
        ...,
        "--out-raw",
        **OUTPUT_FILE_OPTS,
        help="Per-peak output CSV with family + Infernal annotations.",
    ),
    family_csv: Path = typer.Option(
        ...,
        "--family-csv",
        **INPUT_FILE_OPTS,
        help="Lookup CSV (e.g. virus_stats.csv) with at least 'accession' and 'family' columns.",
    ),
    min_treated_sf: float = typer.Option(
        0.0,
        "--min-treated-sf",
        help="Drop peaks with treated_sf below this threshold (applied after "
        "peak_pass==True). Default 0 = no filtering on stall fraction. "
        "treated_sf is k_window / N_treated, the within-tile window stall "
        "fraction (uncorrected for control); use this rather than adj_sf "
        "when you want a filter that doesn't depend on control noise.",
    ),
    min_reads: int = typer.Option(
        0,
        "--min-reads",
        help="Drop peaks with tile_reads_treated below this threshold. Default 0 = no filtering on coverage.",
    ),
    out_collapsed: Optional[Path] = typer.Option(
        None,
        "--out-collapsed",
        **OUTPUT_FILE_OPTS,
        help="Optional per-accession collapsed output CSV (one row per accession, per-peak fields semicolon-joined).",
    ),
    infernal_tbl: List[Path] = typer.Option(
        [],
        "--infernal-tbl",
        help="Infernal .tbl file (repeatable). Adds 4 diagnostic columns per tbl plus any_hit.",
    ),
    evalue: float = typer.Option(0.01, "--evalue", help="Infernal E-value threshold (default 0.001)."),
    downstream: int = typer.Option(
        20,
        "--downstream",
        help="Max nt downstream of stall a hit may start (default 20).",
    ),
    min_downstream_frac: float = typer.Option(
        0.5,
        "--min-downstream-frac",
        help="Min fraction of hit length that must lie downstream of the stall (default 0.5).",
    ),
    additional_columns: str = typer.Option(
        "",
        "--additional-columns",
        help="Comma-separated list of extra peak-CSV columns to carry into the collapsed "
        "output (semicolon-joined per accession alongside the defaults).",
    ),
    min_reps_passing: int = typer.Option(
        0,
        "--min-reps-passing",
        help="Combined-format input only (output of combine_replicates.py): "
        "drop sites where fewer than M replicates had peak_pass=True. "
        "Default 0 = no filter on replicate-pass count (use n_reps_passing column "
        "as informational annotation only).",
    ),
) -> None:
    """Filter, annotate, and (optionally) collapse v7.3 or combine_replicates peaks."""

    df = pd.read_csv(peaks_csv)
    n_input = len(df)

    # ── Detect input format ────────────────────────────────────────────────
    is_combined = "meta_pvalue" in df.columns
    if is_combined:
        print(
            "[postproc] Detected combine_replicates output format "
            "(meta_pvalue column present); filters dispatch to median_treated_sf, "
            "median(tile_reads_treated_*), and n_reps_passing.",
            file=sys.stderr,
        )
        # Required columns for combined-mode filtering.
        for required in ("median_treated_sf", "n_reps_passing", "n_reps_found"):
            if required not in df.columns:
                raise SystemExit(
                    f"[postproc] Combined input is missing '{required}' column — "
                    "is this a current combine_replicates.py output?"
                )
    else:
        for required in ("peak_pass", "treated_sf", "tile_reads_treated"):
            if required not in df.columns:
                raise SystemExit(
                    f"[postproc] Per-replicate input is missing '{required}' column — "
                    "is this a v7.3 filterPeaks output?"
                )

    # ── Filter ─────────────────────────────────────────────────────────────
    if is_combined:
        # Combined-mode equivalents:
        #   peak_pass            ->  n_reps_passing >= min_reps_passing
        #   treated_sf           ->  median_treated_sf
        #   tile_reads_treated   ->  median(tile_reads_treated_R*) across reps
        df_step1 = df[df["n_reps_passing"] >= min_reps_passing]
        n_pass = len(df_step1)
        df_step2 = df_step1[df_step1["median_treated_sf"] >= min_treated_sf]
        n_sf = len(df_step2)
        # Prefer the pre-computed median_tile_reads_treated column if present
        # (current combine_replicates output emits it); otherwise compute on
        # the fly from per-rep tile_reads_treated_* columns (back-compat with
        # older combined CSVs).
        if "median_tile_reads_treated" in df_step2.columns:
            df_step3 = df_step2[df_step2["median_tile_reads_treated"] >= min_reads]
        else:
            rep_read_cols = [c for c in df.columns if c.startswith("tile_reads_treated_")]
            if rep_read_cols:
                median_reads = df_step2[rep_read_cols].median(axis=1, skipna=True)
                df_step3 = df_step2[median_reads >= min_reads]
            else:
                print(
                    "[postproc] WARN: no median_tile_reads_treated or "
                    "tile_reads_treated_* columns found; --min-reads filter "
                    "will not be applied.",
                    file=sys.stderr,
                )
                df_step3 = df_step2
        n_reads = len(df_step3)
        df = df_step3
        print(
            f"[postproc] Filter: {n_input} input -> {n_pass} (n_reps_passing>={min_reps_passing}) "
            f"-> {n_sf} (median_treated_sf>={min_treated_sf}) -> {n_reads} (median tile_reads>={min_reads})",
            file=sys.stderr,
        )
    else:
        df = df[df["peak_pass"] == True]  # noqa: E712 — explicit bool comparison wanted
        n_pass = len(df)
        df = df[df["treated_sf"] >= min_treated_sf]
        n_sf = len(df)
        df = df[df["tile_reads_treated"] >= min_reads]
        n_reads = len(df)
        print(
            f"[postproc] Filter: {n_input} input -> {n_pass} (peak_pass=True) "
            f"-> {n_sf} (treated_sf>={min_treated_sf}) -> {n_reads} (tile_reads>={min_reads})",
            file=sys.stderr,
        )

    # ── Family annotation ───────────────────────────────────────────────────
    df = add_family(df, family_csv)

    # ── Infernal annotation (optional) ──────────────────────────────────────
    if infernal_tbl:
        print(f"[postproc] Annotating Infernal hits from {len(infernal_tbl)} tbl file(s)...", file=sys.stderr)
        df = annotate_infernal(df, infernal_tbl, evalue, downstream, min_downstream_frac)

    # ── Raw output (move taxonomy columns to the front for readability) ─────
    leading = [c for c in TAXONOMY_COLUMNS if c in df.columns]
    raw_cols = leading + [c for c in df.columns if c not in leading]
    df[raw_cols].to_csv(out_raw, index=False)
    print(f"[postproc] Raw output ({len(df)} peaks) -> {out_raw}", file=sys.stderr)

    # ── Collapsed output ────────────────────────────────────────────────────
    if out_collapsed is not None:
        if is_combined:
            # Combined-input default collapse columns.
            default_join = [
                "stall_position",
                "median_treated_sf",
                "mean_treated_sf",
                "meta_pvalue",
                "meta_qvalue",
                "n_reps_found",
                "n_reps_passing",
                "fail_reasons",
            ]
        else:
            default_join = ["stall_position", "adj_sf", "tile_reads_treated", "pvalue", "qvalue"]
        if infernal_tbl:
            # Per-tbl best-E-value columns first (semicolon-joined across the
            # accession's peaks), then the OR-reduced any_hit summary.
            default_join.extend(f"{_tbl_stem(p)}_hit" for p in infernal_tbl)
            default_join.append("any_hit")
        # Default columns may be absent if peaks.csv predates v7.3 — warn but continue.
        missing_default = [c for c in default_join if c not in df.columns]
        if missing_default:
            print(
                f"[postproc] WARN: peaks CSV is missing default collapse columns {missing_default}; "
                "they will be omitted from the collapsed output.",
                file=sys.stderr,
            )
        addl = [c.strip() for c in additional_columns.split(",") if c.strip()]
        missing_addl = [c for c in addl if c not in df.columns]
        if missing_addl:
            print(f"[postproc] WARN: --additional-columns includes unknown column(s): {missing_addl}", file=sys.stderr)
        join_cols = [c for c in default_join if c in df.columns] + [c for c in addl if c in df.columns]
        collapsed = collapse_by_accession(df, join_cols)
        collapsed.to_csv(out_collapsed, index=False)
        print(
            f"[postproc] Collapsed output ({len(collapsed)} accession(s)) -> {out_collapsed}",
            file=sys.stderr,
        )


if __name__ == "__main__":
    app()
