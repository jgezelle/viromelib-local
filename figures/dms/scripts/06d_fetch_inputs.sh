#!/usr/bin/env bash
# 06d_fetch_inputs.sh - run on your LAPTOP. Downloads only the tables 06d_figures.Rmd needs from the
# cluster, keeping the folder layout, so the figures can be drawn locally.
#
# Usage:
#   bash 06d_fetch_inputs.sh                       # defaults below
#   bash 06d_fetch_inputs.sh [user@cluster-host] [local folder]
# Defaults: cluster jgg2144@hpc.c2b2.columbia.edu (the address you ssh to from your Mac),
#           local folder ~/Desktop/viromelib/figures/viromelib_DMS/260903-260930/md2000
#
# Downloads (~300 KB in total):
#   clusters/       pca_scores.tsv pca_loadings.tsv pca_variance.tsv features_zscored.tsv
#   stall_windows/  window_features.tsv window_vs_outside.tsv window_clusters.tsv landscapes.tsv
#                   cm_similarity.tsv taxon_similarity.tsv
#   entropy/ and entropy_test/ (if present)   entropy_windows.tsv entropy_cm_vs_none.tsv
#   plus each folder's run_info.txt
# Not downloaded: work/ (RNAstructure files), entropy_per_nt.tsv, and other large tables.
set -euo pipefail

HOST=${1:-jgg2144@hpc.c2b2.columbia.edu}
VIS=/groups/as6282_gp/scratch_bkup/jgg2144/xrRNASHAPE/DMS/work/06_analysis/virome/260903-260930/md2000
LOCAL=${2:-$HOME/Desktop/viromelib/figures/viromelib_DMS/260903-260930/md2000}
mkdir -p "$LOCAL"

rsync -av --prune-empty-dirs \
  --include='clusters/' --include='stall_windows/' --include='entropy/' --include='entropy_test/' \
  --include='clusters/pca_scores.tsv' --include='clusters/pca_loadings.tsv' \
  --include='clusters/pca_variance.tsv' --include='clusters/features_zscored.tsv' \
  --include='stall_windows/window_features.tsv' --include='stall_windows/window_vs_outside.tsv' \
  --include='stall_windows/window_clusters.tsv' --include='stall_windows/landscapes.tsv' \
  --include='stall_windows/cm_similarity.tsv' --include='stall_windows/taxon_similarity.tsv' \
  --include='entropy*/entropy_windows.tsv' --include='entropy*/entropy_cm_vs_none.tsv' \
  --include='*/run_info.txt' \
  --exclude='*' \
  "$HOST:$VIS/" "$LOCAL/"

echo
echo "Downloaded to $LOCAL. Draw the figures with:"
if [[ -n ${2:-} ]]; then
  echo "  Rscript -e 'rmarkdown::render(\"06d_figures.Rmd\", params = list(vis = \"$LOCAL\"))'"
else
  echo "  Rscript -e 'rmarkdown::render(\"06d_figures.Rmd\")'   (or open it in RStudio and click Knit)"
fi
