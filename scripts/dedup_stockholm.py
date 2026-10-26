#!/usr/bin/env python3
"""
Deduplicate a Stockholm (.stk) alignment file.

LocARNA/mlocarna output frequently contains multiple rows for the same
underlying region (e.g. overlapping sliding-window tiles from the same
peak) that end up as identical aligned sequences once gaps are placed.
This script removes exact duplicate aligned-sequence rows, keeping the
first occurrence of each unique sequence, and rewrites a valid single-block
Stockholm file (updating the #=GF SQ count and preserving #=GC/#=GR lines
such as SS_cons unchanged, since deduplicating rows doesn't change columns).

Usage:
    python3 dedup_stockholm.py <input.stk> <output.stk>

Notes:
- Only handles single-block (non-interleaved) Stockholm files, which is
  what mlocarna produces by default.
- Deduplication is based on the exact aligned sequence string (case-sensitive
  as written, e.g. RNA uppercase U). If you want case-insensitive or
  gap-stripped comparison, see the `key_func` note below.
"""

import sys


def dedup_stockholm(in_path, out_path):
    lines = open(in_path).read().splitlines()

    header = []       # "# STOCKHOLM 1.0"
    gf_lines = []      # #=GF ... lines (comments/metadata)
    seq_lines = []      # (seq_id, seq) in original order
    gc_lines = []      # #=GC / #=GR lines (e.g. SS_cons) — assumed to come after seqs
    footer_seen = False

    for l in lines:
        if l.startswith('# STOCKHOLM'):
            header.append(l)
        elif l.startswith('#=GF'):
            gf_lines.append(l)
        elif l.startswith('#=GC') or l.startswith('#=GR'):
            gc_lines.append(l)
        elif l.strip() == '' or l.strip() == '//':
            continue
        elif l.startswith('#'):
            # any other comment lines, keep as-is
            gf_lines.append(l)
        else:
            parts = l.split(None, 1)
            if len(parts) == 2:
                seq_id, seq = parts
                seq_lines.append((seq_id, seq))

    seen_seqs = set()
    deduped = []
    n_dropped = 0
    for seq_id, seq in seq_lines:
        key = seq  # exact aligned-sequence match; case-sensitive
        if key in seen_seqs:
            n_dropped += 1
            continue
        seen_seqs.add(key)
        deduped.append((seq_id, seq))

    # Update #=GF SQ <n> count if present
    new_gf_lines = []
    for l in gf_lines:
        if l.startswith('#=GF SQ'):
            new_gf_lines.append(f'#=GF SQ {len(deduped)}')
        else:
            new_gf_lines.append(l)

    # Determine column width for alignment (id column) for pretty output
    id_width = max((len(sid) for sid, _ in deduped), default=0)
    gc_tag_lens = []
    for l in gc_lines:
        parts = l.split(None, 2)
        if len(parts) == 3:
            gc_tag_lens.append(len(f'{parts[0]} {parts[1]}'))
    id_width = max(id_width, max(gc_tag_lens, default=0))

    with open(out_path, 'w') as f:
        for l in header:
            f.write(l + '\n')
        for l in new_gf_lines:
            f.write(l + '\n')
        f.write('\n')
        for sid, seq in deduped:
            f.write(f'{sid.ljust(id_width + 2)}{seq}\n')
        for l in gc_lines:
            # Format is "#=GC <feature> <data>" or "#=GR <seqname> <feature> <data>"
            parts = l.split(None, 2)
            if len(parts) == 3:
                mark, feature, data = parts
                tag = f'{mark} {feature}'
                f.write(f'{tag.ljust(id_width + 2)}{data}\n')
            else:
                f.write(l + '\n')
        f.write('//\n')

    print(f"Input rows: {len(seq_lines)}")
    print(f"Duplicate rows removed: {n_dropped}")
    print(f"Unique rows written: {len(deduped)}")
    print(f"Output written to: {out_path}")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.stk> <output.stk>")
        sys.exit(1)
    dedup_stockholm(sys.argv[1], sys.argv[2])
