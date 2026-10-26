#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_igv_session.py — Build an IGV session XML from a set of track files
==========================================================================

Generates a minimal IGV session that loads a genome plus any combination of
5'-end bedGraphs, BAM alignments, GTF feature tracks, and BED feature tracks.
Layout follows the convention used in 3replicates.xml:

  - DataPanel at the top (all bedGraphs stacked together)
  - One Panel per BAM (Coverage + Junctions + Alignment subtracks)
  - FeaturePanel at the bottom (all GTFs and BEDs)

Per-track / per-panel heights and panel-divider fractions are intentionally
omitted on the first pass — IGV fills them in with sensible defaults when
opening a session that lacks them.  If IGV refuses to load the produced
file or the layout looks wrong, this script will need to grow flags for
emitting equal-height Panel attributes; come back to that once we've seen
how IGV handles the minimal session.

All file paths are stored as absolute paths via pathlib's resolve() so the
session is portable to any IGV launch directory on the same machine.

Usage
-----
    python scripts/make_igv_session.py \\
        --genome data/pilot_library_20260319.uniqueIDs.fa \\
        --bedgraph work/11_tracks/control_R1.5p.bedgraph \\
                   work/11_tracks/treated_R1.5p.bedgraph \\
                   work/11_tracks/control_R2.5p.bedgraph \\
                   work/11_tracks/treated_R2.5p.bedgraph \\
                   work/11_tracks/control_R3.5p.bedgraph \\
                   work/11_tracks/treated_R3.5p.bedgraph \\
        --bam work/04_align/fs4/*X1*.bam work/04_align/fs4/*R1*.bam \\
              work/04_align/fs4/*X2*.bam work/04_align/fs4/*R2*.bam \\
              work/04_align/fs4/*X3*.bam work/04_align/fs4/*R3*.bam \\
        --gtf work/06_filter_v7/1_peaks.gtf \\
              work/06_filter_v7/2_peaks.gtf \\
              work/06_filter_v7/3_peaks.gtf \\
        --bed work/05_fragments/pilot_library_20260319.uniqueIDs.bed \\
        --out 3replicates.xml
"""

import argparse
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import List, Optional


# ─────────────────────────────────────────────────────────────────────────────
# Path resolution
# ─────────────────────────────────────────────────────────────────────────────


def resolve(p: Path) -> Path:
    """pathlib equivalent of `realpath`. Errors loudly if the file is missing."""
    rp = p.resolve(strict=False)
    if not rp.exists():
        raise SystemExit(f"[igv-session] File does not exist: {rp}")
    return rp


# ─────────────────────────────────────────────────────────────────────────────
# Track builders.  Attribute keys / clazz strings copied verbatim from the
# reference 3replicates.xml so IGV recognises each track type correctly.
# ─────────────────────────────────────────────────────────────────────────────


def make_bedgraph_track(path: Path, name: str) -> ET.Element:
    return ET.Element(
        "Track",
        {
            "attributeKey": name,
            "clazz": "org.broad.igv.track.DataSourceTrack",
            "id": str(path),
            "name": name,
            "renderer": "BAR_CHART",
            "visible": "true",
            "windowFunction": "mean",
        },
    )


def make_bam_panel(path: Path, panel_name: str) -> ET.Element:
    """One Panel containing Coverage + Junctions + Alignment tracks for a BAM."""
    name = path.name
    panel = ET.Element("Panel", {"name": panel_name})

    ET.SubElement(
        panel,
        "Track",
        {
            "attributeKey": f"{name} Coverage",
            "clazz": "org.broad.igv.sam.CoverageTrack",
            "id": f"{path}_coverage",
            "name": f"{name} Coverage",
            "snpThreshold": "0.2",
            "visible": "true",
        },
    )
    ET.SubElement(
        panel,
        "Track",
        {
            "attributeKey": f"{name} Junctions",
            "clazz": "org.broad.igv.sam.SpliceJunctionTrack",
            "groupByStrand": "false",
            "id": f"{path}_junctions",
            "name": f"{name} Junctions",
            "visible": "false",
        },
    )
    align = ET.SubElement(
        panel,
        "Track",
        {
            "attributeKey": name,
            "clazz": "org.broad.igv.sam.AlignmentTrack",
            "color": "185,185,185",
            "displayMode": "SQUISHED",
            "experimentType": "OTHER",
            "id": str(path),
            "name": name,
            "visible": "true",
        },
    )
    # IGV expects a RenderOptions child on AlignmentTrack to know how to handle dups.
    ET.SubElement(align, "RenderOptions", {"duplicatesOption": "FILTER"})
    return panel


def make_feature_track(path: Path) -> ET.Element:
    name = path.name
    return ET.Element(
        "Track",
        {
            "attributeKey": name,
            "clazz": "org.broad.igv.track.FeatureTrack",
            "groupByStrand": "false",
            "id": str(path),
            "name": name,
            "visible": "true",
        },
    )


# ─────────────────────────────────────────────────────────────────────────────
# Resource type detection
# ─────────────────────────────────────────────────────────────────────────────


_EXT_TO_TYPE = {
    ".bedgraph": "bedgraph",
    ".bg": "bedgraph",
    ".bam": "bam",
    ".gtf": "gtf",
    ".gff": "gff",
    ".gff3": "gff",
    ".bed": "bed",
}


def resource_type(path: Path) -> str:
    """Detect IGV resource type from extension. Lowercase ext only."""
    ext = path.suffix.lower()
    if ext == ".gz":
        ext = path.with_suffix("").suffix.lower()
    return _EXT_TO_TYPE.get(ext, ext.lstrip(".") or "data")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────


def build_session(
    genome: Path,
    bedgraphs: List[Path],
    bams: List[Path],
    gtfs: List[Path],
    beds: List[Path],
    locus: Optional[str],
    version: str,
) -> ET.ElementTree:
    session_attrs = {"genome": str(genome), "version": version}
    if locus:
        session_attrs["locus"] = locus
    session = ET.Element("Session", session_attrs)

    # ── Resources block (flat list of every file referenced) ─────────────────
    resources = ET.SubElement(session, "Resources")
    for p in bedgraphs + bams + gtfs + beds:
        ET.SubElement(resources, "Resource", {"path": str(p), "type": resource_type(p)})

    # ── DataPanel: bedGraphs stacked together at the top ─────────────────────
    if bedgraphs:
        data_panel = ET.SubElement(session, "Panel", {"name": "DataPanel"})
        for p in bedgraphs:
            data_panel.append(make_bedgraph_track(p, name=p.stem))

    # ── One Panel per BAM (Coverage + Junctions + Alignment) ─────────────────
    for i, p in enumerate(bams):
        panel_name = f"Panel_{p.stem}_{i}"
        session.append(make_bam_panel(p, panel_name))

    # ── FeaturePanel: GTFs and BEDs grouped together at the bottom ───────────
    if gtfs or beds:
        feature_panel = ET.SubElement(session, "Panel", {"name": "FeaturePanel"})
        for p in gtfs + beds:
            feature_panel.append(make_feature_track(p))

    # ── HiddenAttributes (matches the reference session) ─────────────────────
    hidden = ET.SubElement(session, "HiddenAttributes")
    for attr_name in ("DATA FILE", "DATA TYPE", "NAME"):
        ET.SubElement(hidden, "Attribute", {"name": attr_name})

    return ET.ElementTree(session)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--genome", required=True, type=Path,
        help="Reference genome FASTA (or .json) for the session.",
    )
    parser.add_argument(
        "--bedgraph", nargs="+", default=[], type=Path,
        help="bedGraph file(s). Stacked into a single DataPanel at the top, in input order.",
    )
    parser.add_argument(
        "--bam", nargs="+", default=[], type=Path,
        help="BAM file(s). Each gets its own Panel (Coverage + Junctions + Alignment), in input order.",
    )
    parser.add_argument(
        "--gtf", nargs="+", default=[], type=Path,
        help="GTF/GFF feature file(s) for the bottom FeaturePanel.",
    )
    parser.add_argument(
        "--bed", nargs="+", default=[], type=Path,
        help="BED feature file(s) for the bottom FeaturePanel.",
    )
    parser.add_argument(
        "--locus", default=None,
        help="Optional default locus, e.g. 'NC_024114.1|2493-2824_inferred_3UTR:0-196'. "
        "If omitted, IGV picks one on first open.",
    )
    parser.add_argument(
        "--version", default="8",
        help="IGV session schema version to write (default 8, matching 3replicates.xml).",
    )
    parser.add_argument(
        "--out", required=True, type=Path,
        help="Output session XML path.",
    )
    args = parser.parse_args()

    if not (args.bedgraph or args.bam or args.gtf or args.bed):
        raise SystemExit("[igv-session] Must specify at least one of --bedgraph, --bam, --gtf, --bed.")

    genome = resolve(args.genome)
    bedgraphs = [resolve(p) for p in args.bedgraph]
    bams = [resolve(p) for p in args.bam]
    gtfs = [resolve(p) for p in args.gtf]
    beds = [resolve(p) for p in args.bed]

    tree = build_session(
        genome=genome,
        bedgraphs=bedgraphs,
        bams=bams,
        gtfs=gtfs,
        beds=beds,
        locus=args.locus,
        version=args.version,
    )

    # Pretty-print so the file is human-readable; matches the reference style.
    ET.indent(tree, space="    ")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    tree.write(str(args.out), encoding="UTF-8", xml_declaration=True)

    print(
        f"[igv-session] {len(bedgraphs)} bedGraph(s), {len(bams)} BAM(s), "
        f"{len(gtfs)} GTF(s), {len(beds)} BED(s)  ->  {args.out}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
