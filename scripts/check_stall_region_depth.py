#!/usr/bin/env python3
"""
Check whether the ACTUAL xrRNA-sized region (default: stall_site to
stall_site + 80 nt, i.e. the structured element itself, not the broader
depth-defined window) has real sequencing depth.

The earlier version of this check just asked "does the stall site fall
somewhere inside the wide depth-window" -- with windows running
200-350+ nt and real xrRNAs only ~70+/-10 nt, that's nearly vacuous: a
wide window containing one point says nothing about whether the actual
structured element has coverage. This version reads each target's own
_profile.txt directly and checks depth across only the biologically
relevant span.

IMPORTANT CAVEAT: this assumes the accession's coordinate numbering
runs in the same direction as the ShapeMapper target's own internal
Nucleotide numbering (i.e., no strand flip happened when the reference
FASTA was extracted). If any of these features were pulled from the
minus strand and reverse-complemented during extraction, position
math here will be silently wrong for those specific targets. Worth
spot-checking a known case (e.g. Zika, where the real xrRNA structure
is published) before trusting this at scale.

Usage:
    python3 check_stall_region_depth.py \
        --master 123_peaks_annotated_with_dms.csv \
        --profiles-glob "$SM_OUT/virome_pooled_Xrn1passing/*_profile.txt" \
        --xrrna-length 80 \
        --min-depth 5000 \
        --out stall_region_depth_check.csv
"""
import argparse
import csv
import glob
import os
import re
import sys
import pandas as pd


def parse_accession_header(accession):
    m = re.match(r'^([^|]+)\|(\d+)-(\d+)', accession)
    if not m:
        return None, None, None
    return m.group(1), int(m.group(2)), int(m.group(3))


def load_profile_depths(path):
    """Returns list of (position, mod_depth, unt_depth), 1-indexed positions matching file order."""
    out = []
    with open(path, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            try:
                out.append((int(row['Nucleotide']),
                            int(row['Modified_read_depth']),
                            int(row['Untreated_read_depth'])))
            except (KeyError, ValueError):
                continue
    return out


def build_profile_index(profiles_glob, run_prefix):
    """Map safe_id (target name minus run prefix) -> profile file path."""
    index = {}
    for path in glob.glob(profiles_glob):
        base = os.path.basename(path).replace('_profile.txt', '')
        safe_id = base.replace(run_prefix, '', 1)
        index[safe_id] = path
    return index


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", required=True)
    ap.add_argument("--profiles-glob", required=True, help="glob for *_profile.txt files")
    ap.add_argument("--run-prefix", default="virome_pooled_Xrn1passing_")
    ap.add_argument("--xrrna-length", type=int, default=80,
                     help="length downstream of the stall site to check (default 80, matching ~70+/-10nt)")
    ap.add_argument("--min-depth", type=int, default=5000)
    ap.add_argument("--min-pct-covered", type=float, default=80.0,
                     help="pct of the xrRNA region that must clear min-depth to count as 'covered'")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.master)
    required = ['accession', 'n_reps_passing', 'stall_site_genomic_position']
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(f"ERROR: master missing columns: {missing}", file=sys.stderr)
        sys.exit(1)

    profile_index = build_profile_index(args.profiles_glob, args.run_prefix)
    print(f"Found {len(profile_index)} profile files matching the glob", file=sys.stderr)

    passing = df[df['n_reps_passing'] == 3].dropna(subset=['stall_site_genomic_position'])

    results = []
    for accession, group in passing.groupby('accession'):
        base_acc, acc_start, acc_end = parse_accession_header(accession)
        if acc_start is None:
            continue

        # find this accession's profile file via the namemap-free safe_id pattern:
        # accession -> safe_id is just '|' replaced with '_'
        safe_id = accession.replace('|', '_')
        profile_path = profile_index.get(safe_id)
        if profile_path is None:
            results.append({
                "accession": accession, "stall_site_genomic": None,
                "xrrna_region_relative": None, "pct_region_covered": None,
                "status": "no_profile_file_found"
            })
            continue

        depths = load_profile_depths(profile_path)
        depth_by_pos = {p: (m, u) for p, m, u in depths}

        for stall_genomic in group['stall_site_genomic_position'].unique():
            stall_genomic = int(stall_genomic)
            stall_rel = stall_genomic - acc_start + 1  # position within this accession's own sequence
            region_start = max(1, stall_rel)
            region_end = stall_rel + args.xrrna_length

            region_positions = [p for p in range(region_start, region_end + 1) if p in depth_by_pos]
            if not region_positions:
                results.append({
                    "accession": accession, "stall_site_genomic": stall_genomic,
                    "xrrna_region_relative": f"{region_start}-{region_end}",
                    "pct_region_covered": 0.0, "status": "region_outside_reference_bounds"
                })
                continue

            n_covered = sum(1 for p in region_positions
                             if depth_by_pos[p][0] >= args.min_depth and depth_by_pos[p][1] >= args.min_depth)
            pct_covered = 100 * n_covered / len(region_positions)
            status = "COVERED" if pct_covered >= args.min_pct_covered else "NOT_COVERED -- xrRNA region itself lacks depth"

            results.append({
                "accession": accession, "stall_site_genomic": stall_genomic,
                "xrrna_region_relative": f"{region_start}-{region_end}",
                "pct_region_covered": round(pct_covered, 1), "status": status
            })

    out_df = pd.DataFrame(results)
    out_df.to_csv(args.out, index=False)
    print(f"\nChecked {len(out_df)} accession/stall-site combinations")
    print(out_df['status'].value_counts())
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
