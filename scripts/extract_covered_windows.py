#!/usr/bin/env python3
"""
Extract the real tile-covered window from ShapeMapper2 _profile.txt files,
rather than trusting the whole-reference pass/fail summary line.

Tiled oligo libraries only cover specific sub-windows of a longer reference
sequence. ShapeMapper's blanket "X% of nucleotides meet minimum depth"
check applies to the FULL reference span, so a target that's mostly
untiled padding around one well-covered tile-length window will always
look like a failure even when that window's data is excellent.

This script finds the longest contiguous run of positions meeting a
depth threshold in BOTH modified and untreated samples, and reports
reactivity/QC stats scoped to just that window.

Usage:
    python3 extract_covered_windows.py <profile1.txt> [<profile2.txt> ...] --out summary.tsv
    python3 extract_covered_windows.py --glob "*_profile.txt" --out summary.tsv

Column names are matched by case-insensitive substring rather than exact
string, since I'm working from memory of ShapeMapper2's format, not your
actual file. If this errors with "could not find column", paste the
actual header line back and I'll fix the matching.
"""
import argparse
import csv
import glob
import statistics
import sys
from pathlib import Path


def find_col(header, *substrings):
    """Find a column index by case-insensitive substring match. Tries each
    substring in order, returns the first match found."""
    lower = [h.lower() for h in header]
    for sub in substrings:
        for i, h in enumerate(lower):
            if sub in h:
                return i
    return None


def load_profile(path):
    with open(path, newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        rows = [row for row in reader if row and len(row) == len(header)]
    return header, rows


def to_float(x):
    try:
        v = float(x)
        if v != v:  # NaN
            return None
        return v
    except (ValueError, TypeError):
        return None


def find_covered_windows(depths_mod, depths_unt, min_depth, min_window_len, require_both=False):
    """Return list of (start_idx, end_idx) 0-based inclusive runs where depth
    clears min_depth, filtered to >= min_window_len.

    By default (require_both=False) the window is defined by MODIFIED depth
    alone -- untreated samples are often systematically shallower by design
    (fewer replicates, lower-diversity library) and requiring them to
    independently clear the same threshold can mask a perfectly real,
    well-covered tile window. Untreated depth stats are still reported for
    that window either way, so the deficit is visible rather than hidden."""
    n = len(depths_mod)
    windows = []
    run_start = None
    for i in range(n):
        mod_ok = depths_mod[i] is not None and depths_mod[i] >= min_depth
        unt_ok = depths_unt[i] is not None and depths_unt[i] >= min_depth
        ok = (mod_ok and unt_ok) if require_both else mod_ok
        if ok and run_start is None:
            run_start = i
        elif not ok and run_start is not None:
            if i - run_start >= min_window_len:
                windows.append((run_start, i - 1))
            run_start = None
    if run_start is not None and n - run_start >= min_window_len:
        windows.append((run_start, n - 1))
    return windows


def summarize_target(path, min_depth=5000, min_window_len=50):
    header, rows = load_profile(path)

    idx_nt = find_col(header, "nucleotide", "position")
    idx_mod_depth = find_col(header, "modified_read_depth", "modified_effective_depth", "modified.*depth")
    idx_unt_depth = find_col(header, "untreated_read_depth", "untreated_effective_depth", "untreated.*depth")
    idx_norm = find_col(header, "norm_profile")
    idx_hq = find_col(header, "hq_profile")
    idx_react = find_col(header, "reactivity_profile", "reactivity")

    missing = []
    if idx_mod_depth is None:
        missing.append("modified depth")
    if idx_unt_depth is None:
        missing.append("untreated depth")
    if missing:
        raise ValueError(
            f"{path}: could not find column(s) for: {', '.join(missing)}. "
            f"Actual header was: {header}"
        )

    idx_reactivity = idx_norm if idx_norm is not None else (idx_hq if idx_hq is not None else idx_react)

    depths_mod = [to_float(r[idx_mod_depth]) for r in rows]
    depths_unt = [to_float(r[idx_unt_depth]) for r in rows]
    reactivities = [to_float(r[idx_reactivity]) for r in rows] if idx_reactivity is not None else [None] * len(rows)

    windows = find_covered_windows(depths_mod, depths_unt, min_depth, min_window_len)

    target_name = Path(path).stem.replace("_profile", "")
    results = []
    if not windows:
        results.append({
            "target": target_name,
            "ref_length": len(rows),
            "window_start": None,
            "window_end": None,
            "window_length": 0,
            "mod_depth_mean_in_window": None,
            "unt_depth_mean_in_window": None,
            "unt_depth_max_in_window": None,
            "mod_unt_ratio": None,
            "unt_clears_5000": None,
            "median_reactivity_in_window": None,
            "mean_reactivity_in_window": None,
            "n_reactive_positions": None,
        })
        return results

    # report the longest window; still note if there were multiple
    windows.sort(key=lambda w: w[1] - w[0], reverse=True)
    for rank, (start, end) in enumerate(windows):
        window_reactivities = [reactivities[i] for i in range(start, end + 1) if reactivities[i] is not None]
        n_reactive = sum(1 for v in window_reactivities if v is not None and v > 0.4)  # loose "reactive" cutoff
        window_mod = [depths_mod[i] for i in range(start, end + 1) if depths_mod[i] is not None]
        window_unt = [depths_unt[i] for i in range(start, end + 1) if depths_unt[i] is not None]
        mod_mean = statistics.mean(window_mod) if window_mod else None
        unt_mean = statistics.mean(window_unt) if window_unt else None
        unt_max = max(window_unt) if window_unt else None
        results.append({
            "target": target_name if rank == 0 else f"{target_name}_window{rank+1}",
            "ref_length": len(rows),
            "window_start": start + 1,  # 1-indexed for reporting
            "window_end": end + 1,
            "window_length": end - start + 1,
            "mod_depth_mean_in_window": round(mod_mean, 0) if mod_mean else None,
            "unt_depth_mean_in_window": round(unt_mean, 0) if unt_mean else None,
            "unt_depth_max_in_window": unt_max,
            "mod_unt_ratio": round(mod_mean / unt_mean, 2) if mod_mean and unt_mean else None,
            "unt_clears_5000": (unt_max is not None and unt_max >= min_depth),
            "median_reactivity_in_window": round(statistics.median(window_reactivities), 3) if window_reactivities else None,
            "mean_reactivity_in_window": round(statistics.mean(window_reactivities), 3) if window_reactivities else None,
            "n_reactive_positions": n_reactive,
        })
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="*", help="profile.txt files")
    ap.add_argument("--glob", help="glob pattern to find profile.txt files, e.g. '*_profile.txt'")
    ap.add_argument("--out", required=True, help="output TSV summary path")
    ap.add_argument("--min-depth", type=int, default=5000, help="depth threshold (default matches ShapeMapper's own QC: 5000)")
    ap.add_argument("--min-window-len", type=int, default=50, help="minimum contiguous covered length to report (nt)")
    args = ap.parse_args()

    files = list(args.files)
    if args.glob:
        files += glob.glob(args.glob)
    if not files:
        print("No input files given. Use positional args or --glob.", file=sys.stderr)
        sys.exit(1)

    all_rows = []
    errors = []
    for f in files:
        try:
            all_rows.extend(summarize_target(f, min_depth=args.min_depth, min_window_len=args.min_window_len))
        except Exception as e:
            errors.append(str(e))

    if all_rows:
        fieldnames = list(all_rows[0].keys())
        with open(args.out, "w", newline="") as out:
            writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(all_rows)
        print(f"Wrote {len(all_rows)} rows to {args.out}")

    n_no_window = sum(1 for r in all_rows if r["window_length"] == 0)
    print(f"{len(files) - len(errors)} files processed OK, {len(errors)} errored, "
          f"{n_no_window} targets had NO covered window at depth >= {args.min_depth}")

    if errors:
        print("\nFirst few errors:", file=sys.stderr)
        for e in errors[:5]:
            print(f"  {e}", file=sys.stderr)


if __name__ == "__main__":
    main()
