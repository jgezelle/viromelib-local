#!/usr/bin/env bash
# Put BOTH laptop infernal trees into one defined layout, organised by DATASET OF ORIGIN.
# Supersedes migrate_cm_organization_261027.sh, which grew stage by stage and moved flaviLib
# builds into the viromelib tree.
#
# ---------------------------------------------------------------------------------------------
# TARGET LAYOUT
#
#   ~/Desktop/viromelib/work/infernal/          SHARED model library + VIROMELIB-dataset builds
#     cm/                 every .cm, one copy each, with MANIFEST.tsv
#       cm_set_<yymmdd>/    the dated set used to annotate a run (git-tracked)
#       seeds/              seed alignments for PROVIDED / shared models (1a 1b 2 3 dumbbell)
#     build/<model>/      viromelib-derived builds: class4, picorna, cl5
#       seed/ iter_vN/ rscape/ tblout/ NOTES.md
#     fa/                 FASTAs, including RNAvirus_CompleteAssem_No_CoV-2.fa + .ssi
#     archive/            retired work, kept for provenance
#
#   ~/Desktop/viromelib/flaviLib/work/infernal/ FLAVILIB-dataset builds only
#     build/<model>/      cand_FJTB_261026, cand_PLLG_260915, class1b, class1c
#       seed/ iter_vN/ rscape/ tblout/ NOTES.md
#     fa/                 per-genus pools, the TRUE-ordered library reference
#     archive/            Sept working trees and retired annotation outputs
#
#   Cluster $ROOT/work/06_infernal/             APPLICATION only -- NOT touched by this script
#
# THREE RULES
#   1. A .cm lives in the shared cm/ -- the same finished models annotate both datasets, and two
#      copies is what produced the Class_1a vs class1a NAME divergence found 2026-10-27.
#   2. A BUILD lives with the dataset whose stalls produced it.
#   3. A seed alignment for a model that was PROVIDED (not built here) lives in cm/seeds/.
# ---------------------------------------------------------------------------------------------
#
# CHECK BY DEFAULT. Nothing moves until --apply. Nothing is overwritten or deleted: a move whose
# destination exists is reported and skipped. One exception, and only with --apply:
# cm/_from_cluster_261027 is removed, because it was a staging copy verified byte-identical to
# cm/cm_set_261021 on 2026-10-27.
#
# Calibration note: cmcalibrate writes its E-value parameters INTO the .cm (the ECM* lines), so a
# model is self-contained and renaming or moving it loses nothing. .i1* files are cmpress output
# for cmscan and are not needed by cmsearch; none exist on either machine.
#
# Usage
#   bash infernal_layout_261027.sh                  # report what would change
#   bash infernal_layout_261027.sh --apply
#   bash infernal_layout_261027.sh --stage <name>   # one stage only
#
# Stages, in the order --apply runs them:
#   bydataset   flaviLib builds out of the viromelib tree
#   models      every .cm into the shared cm/; seeds of provided models into cm/seeds/
#   perbuild    per-model stockholm/ rscape/ tblout/ artifacts into their build/<model>/
#   cmalign     cmalign/ -> build/class4/iter_vN/ (it is entirely the class4 chain)
#   archive     retired annotation outputs and superseded runs into archive/
#   empties     remove directories the moves left empty
#   report      what is left unplaced, and what to write by hand
set -u

VIR="$HOME/Desktop/viromelib/work/infernal"
FLA="$HOME/Desktop/viromelib/flaviLib/work/infernal"
APPLY=0
ONLY=""

while [ $# -gt 0 ]; do
  case "$1" in
    --apply) APPLY=1; shift;;
    --stage) ONLY="$2"; shift 2;;
    --vir) VIR="$2"; shift 2;;
    --fla) FLA="$2"; shift 2;;
    *) echo "unknown arg: $1" >&2; exit 1;;
  esac
done
[ -d "$VIR" ] || { echo "no such directory: $VIR" >&2; exit 1; }
[ -d "$FLA" ] || { echo "no such directory: $FLA" >&2; exit 1; }

want() { [ -z "$ONLY" ] || [ "$ONLY" = "$1" ]; }
say()  { printf '%s\n' "$*"; }
short() { printf '%s' "${1/#$HOME/~}"; }

# Silent in check mode: `move` prints the readable src -> dst summary, so echoing the raw
# mkdir/mv with absolute paths as well made the report unusable.
run() {
  if [ "$APPLY" = 1 ]; then "$@" || say "      FAILED: $*"; fi
}

# move SRC DST -- never clobber, never fail the script
move() {
  local src="$1" dst="$2"
  [ -e "$src" ] || return 0
  if [ -e "$dst" ]; then say "      SKIP, destination exists: $(short "$dst")"; return 0; fi
  run mkdir -p "$(dirname "$dst")"
  run mv "$src" "$dst"
  say "      $(short "$src")"
  say "        -> $(short "$dst")"
}

# move_glob DST_DIR PATTERN...  -- each match keeps its basename
move_glob() {
  local dst="$1"; shift
  local f
  for f in "$@"; do
    [ -e "$f" ] || continue
    move "$f" "$dst/$(basename "$f")"
  done
}

# ============================================================ 1. builds by dataset of origin
if want bydataset; then
say ""
say "=== builds by dataset of origin ==="
say ""
say "  flaviLib-derived -> $(short "$FLA")/build/"
for b in cand_FJTB_261026 cand_PLLG_260915 class1b class1c; do
  move "$VIR/build/$b" "$FLA/build/$b"
done
say ""
say "  staying in $(short "$VIR")/build/ : class4, picorna, cl5 (viromelib-derived)"
say "  class4 has zero hits across all 439 flaviLib accessions, so it is viromelib's."
fi

# ============================================================ 2. the shared model library
if want models; then
say ""
say "=== shared model library -> $(short "$VIR")/cm/ ==="
say ""
say "  candidate models out of their build directories:"
move "$FLA/build/cand_FJTB_261026/cand_FJTB_261026.cm" "$VIR/cm/cand_FJTB_261026.cm"
say ""
say "  a .stk does not belong in cm/:"
move "$VIR/cm/class4_v9_20260824.stk" "$VIR/build/class4/seed/class4_v9_20260824.stk"
say ""
say "  seed alignments for PROVIDED models -> cm/seeds/"
move_glob "$VIR/cm/seeds" \
  "$VIR/stockholm/CL_Class_1a_xrRNA_final.sto" \
  "$VIR/stockholm/CL_Class_1b_xrRNA_final.sto" \
  "$VIR/stockholm/CL_Class_2_xrRNA_final.sto" \
  "$VIR/stockholm/cl1a_manual_jgg.stk" \
  "$VIR/stockholm/cl1b_manual_jgg.stk" \
  "$VIR/stockholm/cl2_manual_jgg.stk" \
  "$VIR/stockholm/cl3_manual_jgg.stk" \
  "$VIR/stockholm/cl3_Gezelle_etal_NAR_merge_v0.align.mpkadded.edited.confident.dedup.IDcurated.SS_cons.sto" \
  "$VIR/stockholm/dumbbell_RF00525.sto"
say ""
say "  the staging copy of the cluster's models, verified byte-identical 2026-10-27:"
if [ -d "$VIR/cm/_from_cluster_261027" ]; then
  if [ "$APPLY" = 1 ]; then
    rm -rf "$VIR/cm/_from_cluster_261027" && say "      removed $(short "$VIR")/cm/_from_cluster_261027"
  else
    say "      would REMOVE $(short "$VIR")/cm/_from_cluster_261027  (the only delete in this script)"
  fi
fi
fi

# ============================================================ 3. per-model build artifacts
if want perbuild; then
say ""
say "=== per-model artifacts -> build/<model>/ ==="
say ""
say "  class4: rscape runs, seeds, per-iteration searches"
move "$VIR/rscape/cl4_v7" "$VIR/build/class4/rscape/cl4_v7"
move "$VIR/rscape/cl4_v9" "$VIR/build/class4/rscape/cl4_v9"
move_glob "$VIR/build/class4/seed" \
  "$VIR/stockholm/cl4_manual_jgg.stk" \
  "$VIR/stockholm/cl4_20260720_result.stk" \
  "$VIR/stockholm/20260720_cl4_input_mlocarna.out" \
  "$VIR/fa/20260720_cl4_input_mlocarna.fa" \
  "$VIR/fa/cl4_20260720_aligned.fasta"
for t in "$VIR"/tblout/cl4_v*.tbl "$VIR"/tblout/cl4_20260720*.tbl; do
  move_glob "$VIR/build/class4/tblout" "$t"
done
say ""
say "  picorna / cl5 (2025 exploration, viromelib)"
move_glob "$VIR/build/picorna/seed" \
  "$VIR"/stockholm/picorna*.stk "$VIR"/stockholm/picorna*.sto \
  "$VIR"/stockholm/Picornaviridae.peak_windows.merge.dedup.out.stk
move_glob "$VIR/build/picorna/tblout" "$VIR"/tblout/picorna*.tbl
move_glob "$VIR/build/cl5/seed"       "$VIR/stockholm/cl5_manual_jgg.stk"
move_glob "$VIR/build/picorna/seed"   "$VIR/coords/picornahit1_coords.txt"
move_glob "$VIR/build/cl5/seed"       "$VIR/coords/cl5coords.txt"
fi

# ============================================================ 3b. cmalign/ is the class4 chain
# From 06_infernal.md: every file in cmalign/ is one step of the class4 v1->v9 iteration. The
# vN_sig_* sets are the significant hits FROM vN, used to build vN+1; they keep their own version
# label so the filename still says which model produced them.
if want cmalign; then
say ""
say "=== cmalign/ -> build/class4/iter_vN/  (the class4 v1->v9 chain) ==="
say ""
say "  v3: the first cmalign of cl4_20260720-2 hits"
move_glob "$VIR/build/class4/iter_v3" \
  "$VIR/cmalign/hits.fetchlist" "$VIR/cmalign/hits.fa" "$VIR/cmalign/hits_aligned.stk"
say ""
say "  v4: diverse known members added by efetch, then merged with the v3 hits"
move_glob "$VIR/build/class4/iter_v4" \
  "$VIR/cmalign/new_accessions.txt" "$VIR/cmalign/new_seqs.fa" "$VIR/cmalign/new_seqs.tbl" \
  "$VIR/cmalign/new_hits.fetchlist" "$VIR/cmalign/new_hits.fa" \
  "$VIR/cmalign/combined_hits.fa" "$VIR/cmalign/combined_aligned.stk"
say ""
say "  v4..v8: each iteration's significant hits, realigned to build the next model"
for v in 4 5 6 7 8; do
  move_glob "$VIR/build/class4/iter_v$v" \
    "$VIR/cmalign/v${v}_sig_hits.fetchlist" \
    "$VIR/cmalign/v${v}_sig_hits.fa" \
    "$VIR/cmalign/v${v}_sig_aligned.stk"
done
say ""
if [ "$APPLY" = 1 ]; then
  say "  left in cmalign/ -- not named in 06_infernal.md, place these by hand:"
else
  say "  cmalign/ currently holds (nothing has moved yet in check mode; after --apply only"
  say "  files not named in 06_infernal.md remain):"
fi
if [ -d "$VIR/cmalign" ]; then ls -A "$VIR/cmalign" 2>/dev/null | sed 's/^/      /'; fi
fi

# ============================================================ 4. archive
if want archive; then
say ""
say "=== archive ==="
say ""
say "  retired annotation grouping -- superseded by the peaks_annotated per-model columns"
move "$VIR/groups_out"   "$VIR/archive/groups_out"
move "$VIR/groups_out_R" "$VIR/archive/groups_out_R"
move "$FLA/groups"       "$FLA/archive/groups"
say ""
say "  csv/ is the INPUT side of that same retired pipeline (infernal.md):"
say "  buildInfernalRunMap.py writes csv/<date>_infernal_runs.csv, one row per tblout, which"
say "  sortclasses.py consumes to produce the window_hits_long / window_counts_by_class /"
say "  summary CSVs now in archive/groups_out/. joinPASScsv_Infernal.py merged the result with"
say "  the PASS table. All of it is replaced by the per-model columns in peaks_annotated."
move "$VIR/csv" "$VIR/archive/csv_runmaps"
say ""
say "  superseded runs"
move "$VIR/260602" "$VIR/archive/260602"
say ""
say "  searches against the 2025 runs, re-processed since. The run IDs appear in filenames as"
say "  r366all / r368 / r0371 / r0373, not consistently zero-padded, so match each form."
for t in "$VIR"/tblout/*r366*.tbl "$VIR"/tblout/*r368*.tbl \
         "$VIR"/tblout/*r0371*.tbl "$VIR"/tblout/*r0373*.tbl \
         "$VIR"/tblout/*r371*.tbl "$VIR"/tblout/*r373*.tbl; do
  move_glob "$VIR/archive/tblout_2025runs" "$t"
done
move_glob "$VIR/archive/tblout_2025runs" "$VIR"/tblout/CL_Class_*.cm.tbl
say ""
say "  leftover coordinate lists"
move_glob "$VIR/archive" "$VIR/coords/coords.txt"
say ""
say "  Sept working trees -- PROVENANCE, cited by the master document's path table"
move "$FLA/flaviLib_20260914" "$FLA/archive/flaviLib_20260914"
move "$FLA/flaviLib_20260917" "$FLA/archive/flaviLib_20260917"
fi

# ============================================================ 5. empties
if want empties; then
say ""
say "=== remove directories left empty ==="
for d in "$FLA/cm/PLLG_260915" "$FLA/cm" "$VIR/rscape" "$VIR/coords" \
         "$VIR/stockholm" "$VIR/tblout" "$FLA/stockholm" "$FLA/tblout"; do
  if [ -d "$d" ] && [ -z "$(ls -A "$d" 2>/dev/null)" ]; then
    run rmdir "$d"
    say "      remove empty $(short "$d")"
  elif [ -d "$d" ]; then
    say "      NOT empty, left alone: $(short "$d")  ($(ls -A "$d" | wc -l | tr -d ' ') entries)"
  fi
done
fi

# ============================================================ 6. report
if want report; then
say ""
say "=== what is left, and what to write by hand ==="
say ""
for t in "$VIR" "$FLA"; do
  say "  $(short "$t")"
  for d in "$t"/*/; do
    [ -d "$d" ] || continue
    say "      $(basename "$d")/   $(ls -A "$d" 2>/dev/null | wc -l | tr -d ' ') entries"
  done
  say ""
done
say "  a. cm/MANIFEST.tsv does not exist. Build it with NAME as a column, since cmsearch writes"
say "     NAME into the tblout and that is where the only real divergence was found:"
say "       cd $(short "$VIR")/cm"
say "       for f in *.cm; do printf '%s\\t%s\\t%s\\n' \"\$f\" \"\$(grep -m1 '^NAME' \$f | awk '{print \$2}')\" \"\$(md5 -q \$f)\"; done"
say ""
say "  b. NOTES.md in each build/<model>/ -- why each iteration was made, hit counts."
say "     Content for cand_FJTB_261026 and cand_PLLG_260915 is in the master document's run log."
say ""
say "  c. Update the master document's path table for anything this script moved."
say ""
say "  d. LAYOUT.txt describes both laptop infernal trees as CM DISCOVERY. Replace with the"
say "     three rules at the top of this script."
say ""
say "  e. git status in BOTH repos. work/infernal/* is git-ignored, and the tracked set is only"
say "     cm/cm_set_261021/*, the candidate models, and each build's seed/ and final .stk."
fi

say ""
if [ "$APPLY" = 1 ]; then say "APPLIED."; else say "CHECK ONLY. Re-run with --apply."; fi
