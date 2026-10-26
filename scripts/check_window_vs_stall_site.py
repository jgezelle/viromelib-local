#!/usr/bin/env python3
"""
Check whether each DMS-covered window (found by depth alone) actually
overlaps the real Xrn1 stall site that defined the accession as a passing
candidate in the first place.

A DMS window can land on sequencing depth for the WRONG tile within a
merged accession span -- this script catches that by comparing:

    accession header's own coordinates (e.g. "NC_015843.2|10322-10990")
      -> gives the genomic offset to convert window_start/end into
         genomic coordinates
    vs.
    stall_site_genomic_position (from the master spreadsheet, for the
      peak row(s) with n_reps_passing == 3 that made this accession pass)

If the stall site does NOT fall inside the DMS window (+/- a small
buffer), the reactivity data for that target is likely NOT measuring
the actual candidate structure, and should be flagged rather than
trusted at face value.

Usage:
    python3 check_window_vs_stall_site.py \
        --master 123_peaks_annotated_with_dms.csv \
        --buffer 10 \
        --out window_stall_overlap_check.csv
"""
import argparse
import re
import sys
import pandas as pd


def parse_accession_header(accession):
    """'NC_015843.2|10322-10990_inferred_3UTR' -> (base_acc, start, end)"""
    m = re.match(r'^([^|]+)\|(\d+)-(\d+)', accession)
    if not m:
        return None, None, None
    return m.group(1), int(m.group(2)), int(m.group(3))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", required=True, help="output of merge_dms_into_master.py (has dms_ columns)")
    ap.add_argument("--buffer", type=int, default=10,
                     help="allow the stall site to be within this many nt of the window edge and still count as overlapping")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.master)

    required = ['accession', 'dms_status', 'dms_window_start', 'dms_window_end',
                'n_reps_passing', 'stall_site_genomic_position']
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(f"ERROR: master is missing expected column(s): {missing}", file=sys.stderr)
        print(f"Available columns: {list(df.columns)}", file=sys.stderr)
        sys.exit(1)

    # only check targets that actually have a DMS window
    has_window = df[df['dms_status'].isin(['good', 'marginal'])].copy()

    results = []
    for accession, group in has_window.groupby('accession'):
        base_acc, acc_start, acc_end = parse_accession_header(accession)
        if acc_start is None:
            continue

        row0 = group.iloc[0]
        wstart, wend = row0['dms_window_start'], row0['dms_window_end']
        if pd.isna(wstart) or pd.isna(wend):
            continue
        wstart, wend = int(wstart), int(wend)

        # convert window (relative to the accession's own sequence, 1-indexed)
        # to genomic coordinates using the accession's own header start
        genomic_window_start = acc_start + wstart - 1 - args.buffer
        genomic_window_end = acc_start + wend - 1 + args.buffer

        # only the peaks that made this accession PASS (n_reps_passing == 3)
        passing_peaks = group[group['n_reps_passing'] == 3]
        stall_sites = passing_peaks['stall_site_genomic_position'].dropna().unique()

        if len(stall_sites) == 0:
            overlap_status = "no_passing_stall_site_found"
        else:
            overlaps = [genomic_window_start <= s <= genomic_window_end for s in stall_sites]
            overlap_status = "OVERLAPS" if any(overlaps) else "NO_OVERLAP -- window may be wrong tile"

        results.append({
            "accession": accession,
            "dms_status": row0['dms_status'],
            "window_relative": f"{wstart}-{wend}",
            "window_genomic": f"{genomic_window_start}-{genomic_window_end}",
            "stall_sites_genomic": ';'.join(str(int(s)) for s in stall_sites),
            "overlap_status": overlap_status,
        })

    out_df = pd.DataFrame(results)
    out_df.to_csv(args.out, index=False)

    print(f"Checked {len(out_df)} accessions with a DMS window")
    print(out_df['overlap_status'].value_counts())
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
