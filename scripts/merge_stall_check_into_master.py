#!/usr/bin/env python3
"""
Merge the stall-region depth check (check_stall_region_depth.py output)
into the master sheet, alongside the existing dms_status columns from
merge_dms_into_master.py.

An accession can have multiple passing peaks/stall sites (e.g. Zika has
two). This collapses to one verdict per accession: COVERED if ANY of its
stall-site regions show real depth, since a single well-covered peak is
still meaningful evidence for that accession.

Adds:
    dms_stall_region_status    "COVERED" / "NOT_COVERED" / "not_checked"
    dms_stall_region_best_pct  the highest pct_region_covered across all
                                stall sites for that accession

Usage:
    python3 merge_stall_check_into_master.py \
        --master 123_peaks_annotated_with_dms.csv \
        --stall-check stall_region_depth_check.csv \
        --out 123_peaks_annotated_with_dms_and_stallcheck.csv
"""
import argparse
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", required=True, help=".csv or .xlsx")
    ap.add_argument("--stall-check", required=True, help="output of check_stall_region_depth.py")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if args.master.lower().endswith(('.xlsx', '.xls')):
        master = pd.read_excel(args.master)
    else:
        master = pd.read_csv(args.master)

    stall = pd.read_csv(args.stall_check)

    # collapse to one row per accession: COVERED if any stall site passed
    def collapse(group):
        covered = (group['status'] == 'COVERED').any()
        best_pct = group['pct_region_covered'].max()
        return pd.Series({
            'dms_stall_region_status': 'COVERED' if covered else 'NOT_COVERED',
            'dms_stall_region_best_pct': best_pct,
        })

    stall_collapsed = stall.groupby('accession').apply(collapse, include_groups=False).reset_index()

    merged = master.merge(stall_collapsed, on='accession', how='left')
    merged['dms_stall_region_status'] = merged['dms_stall_region_status'].fillna('not_checked')

    merged.to_csv(args.out, index=False)

    print(f"Wrote {len(merged)} rows to {args.out}")
    print(merged['dms_stall_region_status'].value_counts())


if __name__ == "__main__":
    main()
