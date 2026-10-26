#!/usr/bin/env python3
"""
Convert a ShapeMapper _profile.txt into a RNAstructure-format .shape file
for ShapeKnots, using the whole reference (not just the covered window --
RNAstructure's -999 "no data" flag lets positions outside real coverage
fold unconstrained, while flanking sequence context still informs the fold).

Prefers Norm_profile; falls back to HQ_profile, then raw Reactivity_profile
if the stricter normalized values are all NaN (as we found happens even
in well-covered regions when ShapeMapper's normalization gate isn't
cleared elsewhere in the target).

Also writes a plain RNA FASTA (T->U) of the same target, since ShapeKnots
needs both.

Usage:
    python3 make_shapeknots_inputs.py \
        --profile virome_pooled_Xrn1passing_NC_035889.1_10329-10808_inferred_3UTR_profile.txt \
        --out-prefix zikv_NC_035889.1
"""
import argparse
import csv
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", required=True)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    with open(args.profile, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        rows = list(reader)

    for col in ['Norm_profile', 'HQ_profile', 'Reactivity_profile']:
        vals = [r[col] for r in rows]
        n_real = sum(1 for v in vals if v != 'nan')
        if n_real > 0:
            chosen_col = col
            break
    else:
        print("ERROR: no usable reactivity column found (all nan)", file=sys.stderr)
        sys.exit(1)

    print(f"Using {chosen_col} ({n_real}/{len(rows)} real values)")

    seq_chars = []
    with open(f"{args.out_prefix}.shape", "w") as shape_out:
        for r in rows:
            pos = r['Nucleotide']
            val = r[chosen_col]
            shape_val = '-999' if val == 'nan' else val
            shape_out.write(f"{pos}\t{shape_val}\n")
            seq_chars.append(r['Sequence'].replace('T', 'U').replace('t', 'u'))

    seq = ''.join(seq_chars)
    with open(f"{args.out_prefix}.fa", "w") as fa_out:
        fa_out.write(f">{args.out_prefix}\n")
        for i in range(0, len(seq), 60):
            fa_out.write(seq[i:i+60] + '\n')

    print(f"Wrote {args.out_prefix}.shape and {args.out_prefix}.fa ({len(seq)} nt)")


if __name__ == "__main__":
    main()
