#!/usr/bin/env python3
"""
Merge the DMS-MaP reactivity summary (pooled_windows.tsv) back onto the
main peaks spreadsheet, keyed by `accession` -- the same join pattern
annotate_peaks_with_clusters.py already uses for cluster_id.

Accepts the master spreadsheet as either .xlsx or .csv.

Adds these columns to the master table:
    dms_status              "good" / "marginal" / "no_window"
    dms_window_start/end    tile-window coordinates found in the DMS profile
    dms_window_length
    dms_mod_depth_mean
    dms_unt_depth_mean
    dms_unt_clears_5000
    dms_mod_unt_ratio
    dms_median_reactivity
    dms_mean_reactivity
    dms_n_reactive_positions

Usage:
    python3 merge_dms_into_master.py \
        --master 123_peaks_annotated20260721-apply-cutoffs20260824.xlsx \
        --windows pooled_windows.tsv \
        --namemap 20260903_pilot_library_20260319_Xrn1resistant_passing_nogenomes.namemap.tsv \
        --run-prefix virome_pooled_Xrn1passing_ \
        --out 123_peaks_annotated_with_dms.csv
"""
import argparse
import csv
import sys
import pandas as pd


def load_namemap(path):
    safe_to_orig = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            safe_to_orig[row['safe_id']] = row['original_id']
    return safe_to_orig


def load_windows(path, safe_to_orig, run_prefix):
    """Returns dict: accession (original, pipe-containing) -> summary dict.
    Only keeps the PRIMARY window per target (skips any secondary '_window2' rows)."""
    by_accession = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            target = row['target']
            last_part = target.rsplit('_', 1)[-1]
            if last_part.startswith('window') and last_part[6:].isdigit():
                continue  # secondary window (e.g. "..._window2") -- skip silently, keep only the primary window
            safe_id = target.replace(run_prefix, '', 1)
            orig = safe_to_orig.get(safe_id)
            if orig is None:
                print(f"WARNING: no namemap entry for {safe_id}, skipping", file=sys.stderr)
                continue

            wlen = int(row['window_length']) if row['window_length'] else 0
            if wlen == 0:
                status = "no_window"
            elif row.get('unt_clears_5000') == 'True':
                status = "good"
            else:
                status = "marginal"

            by_accession[orig] = {
                "dms_status": status,
                "dms_window_start": row.get('window_start') or '',
                "dms_window_end": row.get('window_end') or '',
                "dms_window_length": wlen,
                "dms_mod_depth_mean": row.get('mod_depth_mean_in_window') or '',
                "dms_unt_depth_mean": row.get('unt_depth_mean_in_window') or '',
                "dms_unt_clears_5000": row.get('unt_clears_5000') or '',
                "dms_mod_unt_ratio": row.get('mod_unt_ratio') or '',
                "dms_median_reactivity": row.get('median_reactivity_in_window') or '',
                "dms_mean_reactivity": row.get('mean_reactivity_in_window') or '',
                "dms_n_reactive_positions": row.get('n_reactive_positions') or '',
            }
    return by_accession


def load_master(path):
    if path.lower().endswith('.xlsx') or path.lower().endswith('.xls'):
        return pd.read_excel(path)
    return pd.read_csv(path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", required=True, help="main peaks spreadsheet (.xlsx or .csv)")
    ap.add_argument("--windows", required=True, help="pooled_windows.tsv from extract_covered_windows.py")
    ap.add_argument("--namemap", required=True, help="safe_id -> original_id TSV")
    ap.add_argument("--run-prefix", default="virome_pooled_Xrn1passing_",
                     help="prefix to strip from window targets to get the safe_id")
    ap.add_argument("--out", required=True, help="output CSV path")
    args = ap.parse_args()

    safe_to_orig = load_namemap(args.namemap)
    dms_by_accession = load_windows(args.windows, safe_to_orig, args.run_prefix)

    dms_cols = ["dms_status", "dms_window_start", "dms_window_end", "dms_window_length",
                "dms_mod_depth_mean", "dms_unt_depth_mean", "dms_unt_clears_5000",
                "dms_mod_unt_ratio", "dms_median_reactivity", "dms_mean_reactivity",
                "dms_n_reactive_positions"]

    master = load_master(args.master)
    if 'accession' not in master.columns:
        print(f"ERROR: master has no 'accession' column. Columns found: {list(master.columns)}", file=sys.stderr)
        sys.exit(1)

    dms_df = pd.DataFrame.from_dict(dms_by_accession, orient='index')
    dms_df.index.name = 'accession'
    dms_df = dms_df.reset_index()

    merged = master.merge(dms_df, on='accession', how='left')
    merged['dms_status'] = merged['dms_status'].fillna('not_in_dms_run')

    merged.to_csv(args.out, index=False)

    n_matched = merged['dms_status'].ne('not_in_dms_run').sum()
    print(f"Wrote {len(merged)} rows to {args.out}")
    print(f"{n_matched} rows matched a DMS target ({len(dms_by_accession)} unique accessions had DMS data)")
    print(merged['dms_status'].value_counts())


if __name__ == "__main__":
    main()
