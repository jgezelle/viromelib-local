#!/usr/bin/env bash
# Consolidate CM models and CM-build artifacts into the layout the two infernal READMEs already
# describe. LAPTOP ONLY -- the cluster's work/06_infernal/ is APPLICATION and is not touched.
#
# The target layout (viromelib_infernal_README.md, part A/B):
#
#   ~/Desktop/viromelib/work/infernal/          <- CM DISCOVERY + the shared CM library
#     cm/                  every .cm, one copy, git-tracked, with MANIFEST.tsv
#     build/<cm_name>/     one dir per model being built: seed/ iter_vN/ rscape/ tblout/ NOTES.md
#     db/ archive/
#
#   ~/Desktop/viromelib/flaviLib/work/infernal/ <- flaviLib-SPECIFIC ANALYSES ONLY
#     inputs/<run>/ cl1a_pk/ dumbbell/ archive/
#     (flaviLib_infernal_README.md: "CMs aren't kept here ... CM building and refinement happen
#      there too", i.e. in the shared tree)
#
# DRY RUN BY DEFAULT. Nothing moves until --apply. Nothing is ever overwritten or deleted: a move
# whose destination exists is reported and skipped.
#
# Calibration: cmcalibrate writes its E-value parameters INTO the .cm file (the ECM* lines), so a
# calibrated model is self-contained and renaming or moving it loses nothing. The .i1f/.i1i/.i1m/.i1p
# files are cmpress output for cmscan and are NOT needed by cmsearch. Stage 1 therefore checks for
# ECM lines inside each model, not for sidecar files. (Checked 2026-10-27: all 7 cm_set_261021
# models and both PLLG candidates are calibrated; no .i1* files exist on either machine.)
#
# Usage
#   bash migrate_cm_organization_261027.sh                 # inventory + planned moves
#   bash migrate_cm_organization_261027.sh --apply         # do it
#   bash migrate_cm_organization_261027.sh --stage inventory
set -u

VIR_INF="$HOME/Desktop/viromelib/work/infernal"
FLA_INF="$HOME/Desktop/viromelib/flaviLib/work/infernal"
APPLY=0
STAGE=all

while [ $# -gt 0 ]; do
  case "$1" in
    --apply) APPLY=1; shift;;
    --stage) STAGE="$2"; shift 2;;
    --vir-inf) VIR_INF="$2"; shift 2;;
    --fla-inf) FLA_INF="$2"; shift 2;;
    *) echo "unknown arg: $1" >&2; exit 1;;
  esac
done

[ -d "$VIR_INF" ] || { echo "no such directory: $VIR_INF" >&2; exit 1; }
[ -d "$FLA_INF" ] || { echo "no such directory: $FLA_INF" >&2; exit 1; }

md5of() { md5 -q "$1" 2>/dev/null || md5sum "$1" | awk '{print $1}'; }

say() { printf '%s\n' "$*"; }
run() {
  if [ "$APPLY" = 1 ]; then
    "$@" || say "    FAILED: $*"
  else
    say "    would: $*"
  fi
}

# ------------------------------------------------------------------ move one path, never clobber
move() {
  local src="$1" dst="$2"
  if [ ! -e "$src" ]; then say "    skip (no source): $src"; return; fi
  if [ -e "$dst" ]; then say "    SKIP (destination exists): $dst"; return; fi
  run mkdir -p "$(dirname "$dst")"
  run mv "$src" "$dst"
}

# ------------------------------------------------------------------ move a .cm
# Calibration lives inside the file, so nothing travels with it. Any cmpress sidecars that do exist
# are moved alongside, renamed to match, because cmscan requires the basenames to agree.
move_cm() {
  local src="$1" dst="$2"
  local sbase dbase sdir ddir
  if [ ! -e "$src" ]; then say "    skip (no source): $src"; return; fi
  sbase="$(basename "$src" .cm)"; sdir="$(dirname "$src")"
  dbase="$(basename "$dst" .cm)"; ddir="$(dirname "$dst")"
  grep -q '^ECM' "$src" || say "    WARNING: $sbase has no ECM lines -- not calibrated"
  move "$src" "$dst"
  for ext in i1f i1i i1m i1p; do
    [ -e "$sdir/$sbase.cm.$ext" ] && move "$sdir/$sbase.cm.$ext" "$ddir/$dbase.cm.$ext"
  done
  return 0
}

# ================================================================== 1. inventory
if [ "$STAGE" = all ] || [ "$STAGE" = inventory ]; then
say ""
say "=== 1. INVENTORY: every .cm in both laptop trees ==="
say ""
printf '%-34s %-11s %s\n' "md5(8)" "size" "path"
tmp="$(mktemp)"
for f in $(find "$VIR_INF" "$FLA_INF" -name '*.cm' 2>/dev/null | sort); do
  h="$(md5of "$f")"
  printf '%-34s %-11s %s\n' "${h:0:8}" "$(wc -c < "$f" | tr -d ' ')" "${f/#$HOME/~}"
  printf '%s\t%s\n' "$h" "$f" >> "$tmp"
done
say ""
say "--- duplicate content (same md5, more than one path) ---"
cut -f1 "$tmp" | sort | uniq -d | while read -r h; do
  say "  ${h:0:8}"
  grep "^$h" "$tmp" | cut -f2 | sed "s|^$HOME|~|" | sed 's/^/      /'
done
say "  (none listed above means every .cm is unique content)"
rm -f "$tmp"
say ""
say "--- model NAME (what cmsearch writes into tblout) and calibration state ---"
say "    Two files with the same structure but different NAME are the same model to"
say "    cmsearch's scoring and different models to anything keyed on the name."
printf '%-56s %-16s %s\n' "path" "NAME" "calibrated"
for f in $(find "$VIR_INF" "$FLA_INF" -name '*.cm' 2>/dev/null | sort); do
  nm="$(grep -m1 '^NAME' "$f" | awk '{print $2}')"
  cal=$(grep -q '^ECM' "$f" && echo yes || echo "NO")
  printf '%-56s %-16s %s\n' "$(echo "${f/#$HOME/~}" | tail -c 56)" "${nm:-?}" "$cal"
done
fi

# ================================================================== 2. candidate builds
if [ "$STAGE" = all ] || [ "$STAGE" = builds ]; then
say ""
say "=== 2. CANDIDATE BUILDS -> \$VIR_INF/build/<cand>/ ==="

ORTHO="$FLA_INF/flaviLib_20260917/candidate_class1a/orthoflavi"
FJTB="$VIR_INF/build/cand_FJTB_261026"
say ""
say "  cand_FJTB_261026  (tick-borne 3rd site; was orthoflavi 'cluster 2')"
move     "$ORTHO/orthoflavi_element_261026.fa" "$FJTB/seed/seed.fa"
move     "$ORTHO/element_261026/mlocarna"      "$FJTB/iter_v1/mlocarna"
move     "$ORTHO/element_261026/mlocarna.log"  "$FJTB/iter_v1/mlocarna.log"
move     "$ORTHO/cluster2/mlocarna"            "$FJTB/iter_v0_15seq/mlocarna"
move     "$ORTHO/cluster2/mlocarna.log"        "$FJTB/iter_v0_15seq/mlocarna.log"
move     "$ORTHO/cluster2/rscape"              "$FJTB/iter_v0_15seq/rscape"
move     "$ORTHO/orthoflavi_cluster2_best.fa"  "$FJTB/iter_v0_15seq/seed_15.fa"

PLLG_SRC="$FLA_INF/cm/PLLG_260915"
PLLG="$VIR_INF/build/cand_PLLG_260915"
say ""
say "  cand_PLLG_260915  (was clusterA_candidate / clusterA_refined)"
move_cm  "$PLLG_SRC/clusterA_candidate.cm" "$VIR_INF/cm/cand_PLLG_260915.cm"
move_cm  "$PLLG_SRC/clusterA_refined.cm"   "$VIR_INF/cm/cand_PLLG_260915_cmalign.cm"
# The PLLG files are in stockholm/, not the alignments/ path the 260915 document gives.
move     "$FLA_INF/stockholm/PLLG_mlocarna_out"    "$PLLG/iter_v1/mlocarna"
move     "$FLA_INF/stockholm/PLLG_cmalign_out"     "$PLLG/iter_v2/cmalign"
move     "$FLA_INF/tblout/PLLG_260915"             "$PLLG/tblout"
move     "$FLA_INF/stockholm/PLLG_LGF_unclassified_for_mlocarna.fa" \
         "$PLLG/seed/PLLG_LGF_unclassified_60.fa"
move     "$FLA_INF/stockholm/PLLG_LGF_clade_xrRNA_landscape.csv" \
         "$PLLG/seed/PLLG_LGF_clade_xrRNA_landscape.csv"
say "    left in stockholm/: unclassified_genus_with_mifsud_clade.csv (the clade lookup, used"
say "    beyond PLLG) and cluster_summary.tsv (check which run wrote it before moving)"
fi

# ================================================================== 3. shared CM library
if [ "$STAGE" = all ] || [ "$STAGE" = library ]; then
say ""
say "=== 3. SHARED CM LIBRARY -> \$VIR_INF/cm/ ==="
say ""
say "  class1c: the viromelib README already flags this as needing to move"
move_cm "$FLA_INF/flaviLib_20260914/candidate_short1b/pesti_cmalign/class1c20260917/class1c.cm" \
        "$VIR_INF/cm/class1c.cm"
say ""
say "  CL_Class_* in the flaviLib tree are duplicates of the shared library."
say "  Inventory (stage 1) prints the md5s. If they match \$VIR_INF/cm/, remove the"
say "  flaviLib copies BY HAND after checking -- this script does not delete."
for c in CL_Class_1a_xrRNA_final CL_Class_1b_xrRNA_final CL_Class_2_xrRNA_final; do
  if [ -e "$FLA_INF/cm/$c.cm" ] && [ -e "$VIR_INF/cm/$c.cm" ]; then
    a="$(md5of "$FLA_INF/cm/$c.cm")"; b="$(md5of "$VIR_INF/cm/$c.cm")"
    if [ "$a" = "$b" ]; then say "    $c  IDENTICAL -> flaviLib copy is redundant"
    else say "    $c  DIFFERS  (${a:0:8} vs ${b:0:8}) -- do not remove either; resolve first"; fi
  fi
done
fi

# ================================================================== 4. superseded builds
if [ "$STAGE" = all ] || [ "$STAGE" = archive ]; then
say ""
say "=== 4. SUPERSEDED BUILD ARTIFACTS -> \$VIR_INF/build/<name>/ or archive/ ==="
OLD914="$FLA_INF/flaviLib_20260914"
move_cm "$OLD914/candidate_short1b/pesti_cm/pesti_short1b.cm" \
        "$VIR_INF/build/class1c/iter_pesti_short1b/pesti_short1b.cm"
move_cm "$OLD914/cm_compare/new_cm/class1b_v2.cm" \
        "$VIR_INF/build/class1b/iter_v2/class1b_v2.cm"
move_cm "$OLD914/cm_compare/new_cm/class1b_v3.cm" \
        "$VIR_INF/build/class1b/iter_v3/class1b_v3.cm"
for h in hepaci_hits hepaci_hits_v2 hepaci_unclas; do
  move_cm "$OLD914/superseded_by_class1c/$h.cm" "$VIR_INF/build/class1c/superseded/$h.cm"
done
say ""
say "  NOT MOVED, left as provenance (the campaign files cite these paths):"
say "    $FLA_INF/flaviLib_20260917/.../orthoflavi/orthoflavi_cluster_*"
say "    $FLA_INF/fa/*.fasta          the Aug 27 per-genus pools"
say "    $FLA_INF/stockholm/          seed Stockholms for the established CMs"
say "                                  (cl1a/cl1b/cl2/cl3/cl4/cl5_manual_jgg.stk, CL_Class_*.sto)."
say "                                  These are build artifacts and would belong in"
say "                                  build/<name>/seed/ under the convention, but the run logs"
say "                                  record them here and cite these paths. Decide before moving."
say "    $FLA_INF/tblout/             per-model full-reference searches"
fi

# ================================================================== 5. what to write by hand
if [ "$STAGE" = all ] || [ "$STAGE" = next ]; then
say ""
say "=== 5. BY HAND AFTER THIS ==="
say ""
say "  a. \$VIR_INF/cm/MANIFEST.tsv -- a row per .cm: name, version, md5, origin, status"
say "     (in use / superseded / candidate). Regenerate the md5 column with:"
say "       cd \$VIR_INF/cm && md5 *.cm"
say ""
say "  b. NOTES.md in each build/<name>/ -- why each iteration was made, hit counts."
say "     For cand_FJTB_261026 and cand_PLLG_260915 the content is in"
say "     workflows/CM_CANDIDATES.md and the two campaign files."
say ""
say "  c. LAYOUT.txt: the line listing BOTH laptop infernal directories as 'CM DISCOVERY'"
say "     contradicts flaviLib_infernal_README.md, which reserves the flaviLib tree for"
say "     library-specific analyses. Narrow it, or the next candidate lands in the wrong tree."
say ""
say "  d. The three CL_Class_* copies in the flaviLib tree differ from the shared library ONLY"
say "     in their NAME field (Class_1a vs class1a, etc). The shared-library copies are"
say "     byte-identical to the cluster's, so those are the ones that annotated the 261021"
say "     tables. The flaviLib copies annotated nothing -- remove them once recorded."
say ""
say "  e. git status in BOTH repos before committing -- work/ is git-ignored in places,"
say "     and a moved .cm can silently leave the tracked set (LAYOUT.txt records this"
say "     happening before, recovered with 'git checkout -- work/06_infernal')."
fi

say ""
if [ "$APPLY" = 1 ]; then say "APPLIED."; else say "DRY RUN. Re-run with --apply to move."; fi
