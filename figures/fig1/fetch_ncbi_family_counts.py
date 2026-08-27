#!/usr/bin/env python3
"""
fetch_ncbi_family_counts.py
───────────────────────────
Fetches total RefSeq complete genome counts per viral family from NCBI nuccore.
Used as the grey bar denominator in the Fig 1A-style plot.

Usage:
    python fetch_ncbi_family_counts.py virus_stats.csv ncbi_family_counts.csv --api-key YOUR_KEY

Output: ncbi_family_counts.csv  (columns: family, ncbi_total)

Requirements:
    pip install biopython
"""

import csv, sys, time, argparse
from Bio import Entrez

Entrez.email = "jeaninegezelle@gmail.com"   # <── set your email


def fetch_family_count(family: str) -> int:
    """
    Count RefSeq complete genomes for a viral family in nuccore.
    Stripped to the two filters that are reliably supported:
      - refseq[filter]       : RefSeq sequences only
      - complete[completeness]: complete genomes only
    Family name via [Organism] searches the full taxonomy subtree.
    """
    query = f'"{family}"[Organism] AND refseq[filter] AND complete[completeness]'
    try:
        handle = Entrez.esearch(db="nuccore", term=query, retmax=0)
        record = Entrez.read(handle)
        handle.close()
        return int(record["Count"])
    except Exception as e:
        print(f"  WARNING: failed for {family}: {e}", file=sys.stderr)
        return -1


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv_in",  help="virus_stats.csv")
    parser.add_argument("csv_out", help="output CSV: family, ncbi_total")
    parser.add_argument("--api-key", default=None, help="NCBI API key")
    args = parser.parse_args()

    if args.api_key:
        Entrez.api_key = args.api_key
        sleep_s = 0.11
    else:
        sleep_s = 0.34

    families = sorted({
        row["family"].strip()
        for row in csv.DictReader(open(args.csv_in))
        if row.get("family", "").strip() not in ("", "Unannotated")
    })

    print(f"Fetching counts for {len(families)} families from NCBI nuccore...")
    print(f'Query: "FAMILY"[Organism] AND refseq[filter] AND complete[completeness]\n')

    results = []
    for i, fam in enumerate(families, 1):
        count  = fetch_family_count(fam)
        status = f"{count:6d}" if count >= 0 else "  FAIL"
        print(f"  [{i:3d}/{len(families)}]  {fam:35s}  {status}")
        results.append({"family": fam, "ncbi_total": count})
        time.sleep(sleep_s)

    with open(args.csv_out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["family", "ncbi_total"])
        writer.writeheader()
        writer.writerows(results)

    n_fail = sum(1 for r in results if r["ncbi_total"] < 0)
    n_zero = sum(1 for r in results if r["ncbi_total"] == 0)
    print(f"\nSaved → {args.csv_out}")
    if n_fail:
        print(f"WARNING: {n_fail} families failed — re-run or fill in manually")
    if n_zero:
        print(f"NOTE: {n_zero} families returned 0:")
        for r in results:
            if r["ncbi_total"] == 0:
                print(f"  {r['family']}")


if __name__ == "__main__":
    main()
