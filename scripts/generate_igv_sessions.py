#!/usr/bin/env python3
"""
generate_igv_sessions.py
────────────────────────
Reads the combined peaks xlsx, selects the top N passing tiles by max adj_sf
across replicates, then generates one make_igv_session.py call per accession
(one XML per unique accession/feature).

Usage:
    python generate_igv_sessions.py \
        --peaks 123_peaks_annotated-bobsent-20260825.xlsx \
        --out-dir igv_sessions \
        --top-n 20 \
        --min-reps 2

    Produces:
        igv_sessions/run_igv_sessions.sh   — shell script with all calls
        igv_sessions/top20_peaks.csv       — table of selected peaks

All paths in the generated shell script are relative to the project root
(i.e. run the .sh from the same directory you'd run make_igv_session.py).
"""

import argparse
import csv
import os
import openpyxl
from pathlib import Path


# ── Fixed file paths (matching the make_igv_session.py invocation) ───────────
GENOME  = "data/pilot_library_20260319.uniqueIDs.fa"
BED     = "work/05_fragments/pilot_library_20260319.uniqueIDs.bed"

BEDGRAPHS = [
    "work/11_tracks/control_R1.5p.bedgraph",
    "work/11_tracks/treated_R1.5p.bedgraph",
    "work/11_tracks/control_R2.5p.bedgraph",
    "work/11_tracks/treated_R2.5p.bedgraph",
    "work/11_tracks/control_R3.5p.bedgraph",
    "work/11_tracks/treated_R3.5p.bedgraph",
]

BAMS = [
    "work/04_align/fs4/*X1*k1.mq0.bam",
    "work/04_align/fs4/*R1*k1.mq0.bam",
    "work/04_align/fs4/*X2*k1.mq0.bam",
    "work/04_align/fs4/*R2*k1.mq0.bam",
    "work/04_align/fs4/*X3*k1.mq0.bam",
    "work/04_align/fs4/*R3*k1.mq0.bam",
]

GTFS = [
    "work/06_filter_v7/1_peaks.gtf",
    "work/06_filter_v7/2_peaks.gtf",
    "work/06_filter_v7/3_peaks.gtf",
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--peaks",    required=True,
                   help="Combined peaks xlsx (e.g. 123_peaks_annotated-bobsent-20260825.xlsx)")
    p.add_argument("--out-dir",  default="igv_sessions",
                   help="Output directory for XMLs and shell script (default: igv_sessions)")
    p.add_argument("--top-n",    type=int, default=20,
                   help="Number of top peaks to select (default: 20)")
    p.add_argument("--min-reps", type=int, default=2,
                   help="Minimum n_reps_passing to include (default: 2)")
    p.add_argument("--sheet",    default="123_peaks_annotated",
                   help="Sheet name in xlsx (default: 123_peaks_annotated)")
    return p.parse_args()


def feature_length(accession: str) -> int:
    """
    Infer feature length from the accession/contig name.
    Format: NCBI_ACC|START-END_FEATURE
    Length = END - START (0-based coordinate span in the tiled FASTA).
    """
    try:
        coords = accession.split("|")[1].split("_")[0]   # e.g. "2493-2824"
        start, end = coords.split("-")
        return int(end) - int(start)
    except Exception:
        return 0


def safe_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── Load xlsx ─────────────────────────────────────────────────────────────
    print(f"Loading {args.peaks} ...")
    wb = openpyxl.load_workbook(args.peaks, read_only=True)
    ws = wb[args.sheet]
    rows = list(ws.iter_rows(values_only=True))
    headers = rows[0]
    col = {h: i for i, h in enumerate(headers) if h}

    adj_cols = [h for h in headers if h and "adj_sf" in str(h).lower()]
    print(f"adj_sf columns: {adj_cols}")

    # ── Collect passing tiles, score by max adj_sf ────────────────────────────
    candidates = []
    for r in rows[1:]:
        n_pass = r[col["n_reps_passing"]]
        if n_pass is None or int(n_pass) < args.min_reps:
            continue

        accession  = r[col["accession"]]      # contig name in FASTA
        family     = r[col["family"]] or "Unannotated"
        tile_start = r[col["tile_start"]]
        tile_end   = r[col["tile_end"]]
        n_reps     = int(n_pass)

        adj_vals = [safe_float(r[col[c]]) for c in adj_cols]
        adj_vals = [v for v in adj_vals if v is not None]
        max_adj  = max(adj_vals) if adj_vals else 0.0

        candidates.append({
            "max_adj_sf":  max_adj,
            "accession":   accession,
            "family":      family,
            "tile_start":  tile_start,
            "tile_end":    tile_end,
            "n_reps":      n_reps,
        })

    print(f"Total passing tiles (n_reps>={args.min_reps}): {len(candidates)}")

    # ── Deduplicate: one entry per accession (best tile) ─────────────────────
    best = {}
    for c in candidates:
        acc = c["accession"]
        if acc not in best or c["max_adj_sf"] > best[acc]["max_adj_sf"]:
            best[acc] = c

    ranked = sorted(best.values(), key=lambda x: x["max_adj_sf"], reverse=True)
    top    = ranked[:args.top_n]
    print(f"Unique accessions: {len(best)} | Selecting top {len(top)}")

    # ── Write peak table ──────────────────────────────────────────────────────
    csv_path = out_dir / f"top{args.top_n}_peaks.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "rank","max_adj_sf","accession","family",
            "tile_start","tile_end","n_reps","locus","xml_out"
        ])
        writer.writeheader()
        for rank, peak in enumerate(top, 1):
            flen  = feature_length(peak["accession"])
            locus = f"{peak['accession']}:0-{flen}" if flen else peak["accession"]
            # safe filename: replace special chars
            safe  = peak["accession"].replace("|","_").replace("/","-")
            peak["rank"]    = rank
            peak["locus"]   = locus
            peak["xml_out"] = f"{args.out_dir}/{safe}.xml"
            writer.writerow(peak)
    print(f"Peak table -> {csv_path}")

    # ── Write shell script ────────────────────────────────────────────────────
    sh_path = out_dir / "run_igv_sessions.sh"
    bedgraph_str = " \\\n               ".join(BEDGRAPHS)
    bam_str      = " \\\n          ".join(BAMS)
    gtf_str      = " \\\n          ".join(GTFS)

    lines = [
        "#!/usr/bin/env bash",
        "# Auto-generated by generate_igv_sessions.py",
        f"# Top {args.top_n} peaks by max adj_sf (n_reps_passing >= {args.min_reps})",
        f"# Source: {args.peaks}",
        "",
        f"mkdir -p {args.out_dir}",
        "",
    ]

    for peak in top:
        flen  = feature_length(peak["accession"])
        locus = f"{peak['accession']}:0-{flen}" if flen else peak["accession"]
        safe  = peak["accession"].replace("|","_").replace("/","-")
        xml   = f"{args.out_dir}/{safe}.xml"

        lines += [
            f"# Rank {peak['rank']} | adj_sf={peak['max_adj_sf']:.4f}"
            f" | {peak['family']} | {peak['accession']}",
            f"python scripts/make_igv_session.py \\",
            f"    --genome {GENOME} \\",
            f"    --bedgraph {bedgraph_str} \\",
            f"    --bam {bam_str} \\",
            f"    --gtf {gtf_str} \\",
            f"    --bed {BED} \\",
            f'    --locus "{locus}" \\',
            f"    --out {xml}",
            "",
        ]

    with open(sh_path, "w") as f:
        f.write("\n".join(lines))
    sh_path.chmod(0o755)
    print(f"Shell script -> {sh_path}")

    # ── Print summary ─────────────────────────────────────────────────────────
    print(f"\nTop {len(top)} peaks:")
    print(f"{'Rank':>4}  {'adj_sf':>8}  {'Family':25s}  Accession")
    for p in top:
        print(f"{p['rank']:>4}  {p['max_adj_sf']:>8.4f}  {p['family']:25s}  {p['accession']}")

    print(f"\nRun from project root:")
    print(f"  bash {sh_path}")


if __name__ == "__main__":
    main()
