# Replicate filtering pipeline (v7)

End-to-end recipe for going from per-replicate alignment BAMs to
post-processed, cross-replicate-combined Xrn1 stall peak calls, with
downstream structural clustering and consensus 2-D visualisation.

This is the same skeleton and scripts as `replicates_filtering_v7.3.md`
(both call `filterPeaks_normtoCTRL_v7.3.py` with calibrated Negative
Binomial p-value gating), kept as a separate doc mainly for its own
`outdir` (`work/06_filter_v7` vs `work/06_filter_v7.3`). The two are
meant to stay in sync — check both when changing a filter default or
adding a pipeline step.

The pipeline assumes three biological replicates labelled `1`, `2`, `3`
where each replicate has a *treated* (Xrn1-digested) and a *control*
(undigested) library, named in the BAM files with the patterns `X{r}`
(treated) and `R{r}` (control) — e.g. `X1`, `R1`, `X2`, `R2`, …

## 0. Initialise paths

Set the variables once at the top of your shell session so every later
command picks up the same paths. `tiles_file` is the BED of tile
coordinates (each tile is an `accession|start-end_label` interval in
the pilot library FASTA); peak calling is tile-aware so that 5′-stall
positions are reported in tile-local coordinates.

```bash
replicates=("1" "2" "3")
peak_script="scripts/filterPeaks_normtoCTRL_v7.3.py"
calibrate_script="scripts/calibratePeakNull.py"
outdir="work/06_filter_v7"
fragment_dir="work/05_fragments"
tiles_file="work/05_fragments/pilot_library_20260319.uniqueIDs.bed"
mkdir -p "$outdir"
```

## 1. Build fragment files

Each BAM is collapsed to a per-fragment CSV (one row per read pair)
with the 5′-end coordinate and a few QC flags. Downstream steps only
ever see these CSVs — the BAMs themselves aren't touched again.

```bash
mkdir -p "$fragment_dir"
for f in /groups/as6282_gp/scratch_bkup/jgg2144/xrrna_seq/viromelib/work/04_align/fs4/*.bam; do
    echo "$f"
    python3 scripts/make_fragment_file.py "$f" \
        > "$fragment_dir/$(basename "$f" .bam)"
done
```

## 2. Calibrate the Negative Binomial null

The v7.3 peak filter adds a per-peak p-value gate (G7) on top of the
v7.2 stall-fraction gates. To compute meaningful p-values we need a
null distribution for per-position 5′-end counts under "no stall". A
single global Negative Binomial dispersion α̂ is fit from
low-stall-fraction background positions across the whole library, per
replicate. The calibrator also produces Q-Q plots so you can sanity-
check the fit.

Look at the printed α̂ for each replicate — you'll paste them into the
next step. Background defaults: positions with SF ≤ 0.05, window=3,
bootstrap=200.

```bash
for r in "${replicates[@]}"; do
    python3 "$calibrate_script" \
        "$fragment_dir"/*-X${r}*k1.mq0.csv \
        "$tiles_file" \
        "$outdir/${r}_calibrate" \
        | tee "$outdir/${r}_calibrate.out.txt"
done
```

Outputs per replicate:
- `${r}_calibrate.alpha.txt` — point estimate + bootstrap CI for α̂
- `${r}_calibrate.qq.png` — Binomial vs NB Q-Q diagnostic plot
- `${r}_calibrate.out.txt` — full log

## 3. Call peaks with the calibrated null

Plug the three α̂ values into the array below (one per replicate, in
the same order as `replicates`). The peak caller then runs all G1–G6
filters plus the new G7 NB p-value gate, with per-accession (default)
or global BH FDR scope.

Key flags:
- `--null negbinom --nb-alpha $alpha` — use the calibrated NB null
  (default would be Binomial, which over-rejects)
- `--min-sf-treated 0.25` — [G2] minimum treated stall fraction at the
  peak top (default is 0.05; raised here so SF filtering happens once,
  at peak-calling, rather than being re-applied in step 6)
- `--pvalue-window 3` — count reads in a ±3-nt window around each peak
  for the p-value test (default; matches the calibrator's `--window`)
- `--max-qvalue 0.05` and `--fdr-scope global` (defaults) — BH across
  every candidate peak in the library

```bash
alphas=("2.0885" "1.7196" "2.1398")

for i in "${!replicates[@]}"; do
    r="${replicates[$i]}"
    alpha="${alphas[$i]}"
    echo "Processing replicate: $r  (α̂ = $alpha)"
    python3 "$peak_script" \
        "$fragment_dir"/*X${r}-*k1.mq0.csv \
        "$fragment_dir"/*R${r}-*k1.mq0.csv \
        "$tiles_file" \
        "$outdir/${r}_peaks.csv" \
        --out-accession-csv "$outdir/${r}_accessions.csv" \
        --out-gtf "$outdir/${r}_peaks.gtf" \
        --min-sf-treated 0.25 \
        --null negbinom \
        --nb-alpha "$alpha"
done
```

Each replicate produces:
- `${r}_peaks.csv` — one row per peak that passes all gates, with the
  `peak_id` / `peak_position_id` columns used to join across replicates
- `${r}_accessions.csv` — per-tile summary (counts, top peak)
- `${r}_peaks.gtf` — peaks as a GTF for IGV display

## 4. Annotate with Infernal CM hits

Run Jeanine's reference covariance models (CMs) against the full
pilot library FASTA. The resulting `--tblout` tables get joined to the
peak calls in step 6 so you can flag which peaks overlap a known
xrRNA family member. Run this once; the tables aren't replicate-
specific.

```bash
cms=(
    "CL_Class_1a_xrRNA_final.cm"
    "CL_Class_1b_xrRNA_final.cm"
    "CL_Class_2_xrRNA_final.cm"
    "cl3_NAR_2025.align.mpkadded.edited.confident.dedup.IDcurated.SS_cons.cm"
    "cl4.cm"
)
cmdir="data/jeanine_cm"
ref="data/pilot_library_20260319.uniqueIDs.fa"
tblout="work/08_infernal"
mkdir -p "$tblout"

for cm in "${cms[@]}"; do
    echo "Processing CM: $cm"
    cmsearch --tblout "$tblout/$cm.tbl" "$cmdir/$cm" "$ref"
done
```

## 5. Combine replicates with Stouffer's method

Cross-replicate combining is tile-aware: peaks that fall in the same
tile and within `--cluster-distance` (default 3 nt) get merged, and
their per-replicate p-values are combined via sqrt-N-weighted
Stouffer's method. The `meta_qvalue` is BH-adjusted across all merged
peaks. Missing replicates default to *neutral* (p=0.5) rather than
being dropped or treated as failures — this is configurable via
`--missing`.

```bash
replicate_csvs=$(printf -- "$outdir/%s_peaks.csv " "${replicates[@]}")
python3 scripts/combine_replicates.py $replicate_csvs \
    --labels R1 R2 R3 \
    --out "$outdir/123_peaks_combined_with_seqs.csv" \
    --fasta data/pilot_library_20260319.uniqueIDs.fa \
    --upstream 0
```

Output: `123_peaks_combined_with_seqs.csv` — the "csv with
everything" — with `peak_id` / `peak_position_id` identifiers,
per-replicate read counts, `meta_pvalue`, `meta_qvalue`,
`n_reps_found`, `n_reps_passing`, median/mean treated SF columns, and
(since `--fasta` is passed) `downstream_seq` / `downstream_seq_start` /
`downstream_seq_end` / `downstream_stall_at_offset` /
`downstream_seq_note` — the downstream sequence for *every* row,
unfiltered (this script already emits one row per cross-rep site, so
there's no tiebreak/dedup step to disable). `downstream_seq_note` is
`ok`, `missing_accession`, or `empty_window`. Everything downstream
(post-processing, clustering FASTA extraction) reads this one file —
there's no separate sequence-free `123_peaks_combined.csv`.

## 6. Post-process: annotate + collapse

`post_process_peaks.py` auto-detects whether the input is single-
replicate or combined (it looks for the `meta_pvalue` column).

Annotates each peak with:
- family / genus / species from the virus-stats reference
- `any_hit` (boolean) plus per-CM hit flags from the `--infernal-tbl`
  arguments
- collapses to one row per accession via OR over the hit columns

Filter knobs (per-replicate runs use them; the combined run can leave
defaults to retain all peaks):
- `--min-reads` — minimum treated read count at the peak
- `--min-treated-sf` — minimum treated stall fraction; left at default
  here since step 3's `--min-sf-treated 0.25` already gates on SF at
  peak-calling time — don't double-apply a second, looser SF filter
  post-hoc

```bash
annout="work/10_postproc"
mkdir -p "$annout"

infernal_args=$(printf -- "--infernal-tbl $tblout/%s.tbl " "${cms[@]}")

for r in "${replicates[@]}"; do
    echo "Annotating peaks for replicate: $r"
    python3 scripts/post_process_peaks.py "$outdir/${r}_peaks.csv" \
        --out-raw "$annout/${r}_peaks_annotated.csv" \
        --out-collapsed "$annout/${r}_peaks_per_accession.csv" \
        --min-reads 20 \
        --family-csv data/virus_stats.csv \
        $infernal_args
done

# Same thing for the combined-replicate call, no SF/reads filter so
# every cross-replicate-supported peak survives into the annotated CSV
python3 scripts/post_process_peaks.py "$outdir/123_peaks_combined_with_seqs.csv" \
    --out-raw "$annout/123_peaks_annotated.csv" \
    --out-collapsed "$annout/123_peaks_per_accession.csv" \
    --family-csv data/virus_stats.csv \
    $infernal_args
```

`post_process_peaks.py` only reads the columns it needs, so the
`downstream_seq*` columns from step 5 pass straight through into
`123_peaks_annotated.csv` untouched.

## 7. Extract downstream sequences for clustering

`123_peaks_combined_with_seqs.csv` already has a `downstream_seq` on
every row (step 5), so this step is purely about producing the
*filtered, deduped* FASTA that LocARNA + clustering actually consume
— the candidate xrRNA cores worth structurally comparing. Dedup at
sequence level keeps the FASTA tidy and writes a sidecar TSV mapping
each unique sequence to the list of peak IDs that contributed it.

Defaults: `--max-length 80`, `--min-passing-replicates 2`,
`--min-treated-sf 0.5`, `--min-treated-reads 20`. Pass `--upstream 0`
to match step 5 (no 5′-context buffer before the stall position).

```bash
python3 scripts/extract_downstream_seqs.py \
    "$outdir/123_peaks_combined_with_seqs.csv" \
    --fasta data/pilot_library_20260319.uniqueIDs.fa \
    --out   "$outdir/123_peaks_combined_downstream_seqs_unique.fa" \
    --upstream 0
```

`--out-csv` is also available here if you want to see exactly which
rows survived filtering/tiebreak/dedup into the FASTA (`ok`) vs why
each other row didn't (`failed_filters`, `lost_tiebreak`,
`duplicate_sequence`, `missing_accession`, `empty_window`) — see the
script's `--out-csv` help for details. For "give me a sequence on
every row" you already have that from step 5; this flag is for
auditing the *filtered* view specifically.

## 8. IGV tracks (optional but useful for QC)

Convert the fragment CSVs to 5′-end bedgraphs (one per replicate per
condition) plus a tile-boundary BED. Then build an IGV session XML
pointing at all of them. The `--locus` flag opens IGV pre-positioned
on a hand-picked accession — adjust as needed.

```bash
tracks_dir="work/11_tracks"
mkdir -p "$tracks_dir"

for r in "${replicates[@]}"; do
    python3 scripts/fragments_to_tracks.py "$fragment_dir"/*X$r-*k1.mq0.csv \
        --tile-bed "$tiles_file" \
        --out-bedgraph "$tracks_dir/treated_R$r.5p.bedgraph" \
        --out-tile-bed "$tracks_dir/treated_R$r.tiles.bed" \
        --track-name "treated_R$r"
    python3 scripts/fragments_to_tracks.py "$fragment_dir"/*R$r-*k1.mq0.csv \
        --tile-bed "$tiles_file" \
        --out-bedgraph "$tracks_dir/control_R$r.5p.bedgraph" \
        --out-tile-bed "$tracks_dir/control_R$r.tiles.bed" \
        --track-name "control_R$r"
done

python3 scripts/make_igv_session.py \
    --genome data/pilot_library_20260319.uniqueIDs.fa \
    --bedgraph work/11_tracks/control_R1.5p.bedgraph \
               work/11_tracks/treated_R1.5p.bedgraph \
               work/11_tracks/control_R2.5p.bedgraph \
               work/11_tracks/treated_R2.5p.bedgraph \
               work/11_tracks/control_R3.5p.bedgraph \
               work/11_tracks/treated_R3.5p.bedgraph \
    --bam work/04_align/fs4/*X1*k1.mq0.bam work/04_align/fs4/*R1*k1.mq0.bam \
          work/04_align/fs4/*X2*k1.mq0.bam work/04_align/fs4/*R2*k1.mq0.bam \
          work/04_align/fs4/*X3*k1.mq0.bam work/04_align/fs4/*R3*k1.mq0.bam \
    --gtf  work/06_filter_v7/{1,2,3}_peaks.gtf \
    --bed  work/05_fragments/pilot_library_20260319.uniqueIDs.bed \
    --locus "NC_024114.1|2493-2824_inferred_3UTR:0-196" \
    --out   3replicates_generated.xml
```

## 9. Pairwise LocARNA (for similarity matrix → clustering)

We need every-vs-every pairwise structural similarity to drive the R
clustering. Locarna's `--local-output` + `--sequ-local=1` + `--struct-local=1`
runs local alignment per pair, which is the right choice here because
the pairwise step is just measuring similarity — global mode would
penalise length differences that aren't biologically real.

The `fasta_split` step turns the combined FASTA into one file per
sequence so locarna can iterate pairs. Replace with any
`split-multi-fasta` equivalent if you don't have it.

```bash
mkdir -p "$outdir/seq_split"
cd "$outdir/seq_split"

# Splits a multi-FASTA into single-FASTA files named by sequence ID.
# Output files look like NC_038434.1|10194-10465_3UTR|191.fa
fasta_split ../123_peaks_combined_downstream_seqs.fa

mkdir -p locarna_out
for x in *.fa; do
    for y in *.fa; do
        if [[ "$x" < "$y" ]] || [[ "$x" == "$y" ]]; then
            echo "locarna --stockholm \"locarna_out/${x}--${y}.stk\" --score-components --local-output --sequ-local=1 --struct-local=1 \"$x\" \"$y\" > \"locarna_out/${x}--${y}.out\""
        fi
    done
done > locarna_commands.sh

# Run them — parallelise with GNU parallel or xargs as fits your cluster
bash locarna_commands.sh        # or: parallel -j 16 < locarna_commands.sh
cd -
```

After every pair has run, assemble the per-pair scores into a single
N×N similarity matrix. clustering.R expects normalised scores (it
picks its connected-components threshold from the 0..1 range), which
means passing `--norm-output`. Normalisation is
`S_ij / sqrt(S_ii · S_jj) × min(1, len_ij / 80)` — coverage-weighted
so that short local alignments don't get full credit.

```bash
python3 scripts/build_locarna_matrix.py \
    "$outdir/seq_split/locarna_out" \
    "$outdir/seq_split/locarna_matrix.csv" \
    --norm-output
```

The output `locarna_matrix.csv` is the file `clustering.R` loads at
its Stage 1 (the `matrix_csv` variable at the top of the script).

## 10. Clustering in R

Open `scripts/clustering.R` and step through it. The 14 stages do:

1. **Load** the pairwise locarna scores into an N×N similarity matrix
2. **Normalise** via `S_ij / sqrt(S_ii · S_jj)` with coverage weighting
3. **KDE-based mode detection** on the similarity distribution to find
   natural thresholds
4. **Percolation sweep** across thresholds, picking the connected-
   components threshold (~0.50 for this library) where the giant
   component first stabilises
5. **Connected-components clustering** at the chosen threshold; cluster
   IDs annotated with `cm_class` (which xrRNA family each known member
   belongs to)
6. **Inter-component dendrogram** + force-directed graph for visual QC
7. **Per-component heatmaps** to confirm intra-cluster cohesion
8. **FASTA emission** — one file per cluster into
   `work/13_msa_groups_hybrid/` (these are mlocarna's input in step 11)
9. **UMAP + spectral clustering** as sanity checks (negative result for
   modular clustering — documented in the script comments)
10. **Hybrid partition**: knowns at τ=0.50, unknowns at τ≈0.47, both
    feeding into `final_comp_df` with `K###` / `N###` prefixes (Known
    vs Novel)

The output FASTAs in `work/13_msa_groups_hybrid/<cluster>.fa` are what
mlocarna consumes next.

## 10b. Join cluster IDs back onto the annotated peaks CSV

Once clustering.R has emitted per-cluster FASTAs, every peak that
contributed a sequence has a cluster (`K###_n##_known` for clusters of
known xrRNA family members, `N###_n##_novel` for unclassified ones).
This step folds that assignment back into the annotated peaks CSV so
downstream analysis (figures, tables, filters by family) can pivot on
`cluster_id` directly without going through the FASTAs.

The script walks every `work/13_msa_groups_hybrid/<cluster_id>.fa`,
treats each FASTA header as a `peak_position_id`, and joins the
resulting `peak_position_id → cluster_id` map onto the CSV. Peaks
whose position isn't in any cluster — filtered out by
`extract_downstream_seqs.py`, or singletons too small to cluster —
get an empty `cluster_id`.

```bash
python3 scripts/annotate_peaks_with_clusters.py \
    "$annout/123_peaks_annotated.csv" \
    "$annout/123_peaks_annotated_with_clusters.csv" \
    --msa-groups-dir work/13_msa_groups_hybrid
```

If the same canonical sequence ended up in both a K and an N cluster
(hybrid partition can do this for borderline peaks), the script warns
on stderr and keeps the first-encountered assignment — alphabetic
order means K wins over N. Adjust if you want the opposite priority.

## 11. Multiple alignment with mlocarna

Globally aligns the per-cluster FASTAs. Uses `--probabilistic` for
better-quality alignment and `--consensus-structure alifold` to derive
the consensus secondary structure via RNAalifold. The driver script
is idempotent — it skips clusters that already have a Stockholm
result, so you can re-run after editing clustering.R without
re-aligning everything.

```bash
bash scripts/run_mlocarna.sh
# Defaults: work/13_msa_groups_hybrid → work/14_msa, 4 threads
# Override: bash scripts/run_mlocarna.sh <in_dir> <out_dir> <threads>
```

Output per cluster:
- `work/14_msa/<cluster>/results/result.stk` — Stockholm with
  consensus structure
- `work/14_msa/<cluster>.log` — per-cluster mlocarna log

## 12. Per-cluster MSA + consensus-structure figures

Two complementary visualisations:

**12a. R-based 3-panel MSA figure** (arcs + alignment heatmap +
conservation), one PNG/SVG per cluster:

```bash
Rscript scripts/visualize_msa.R
# Outputs land in work/15_msa_plots/<cluster>.{png,svg}
```

**12b. R2R consensus 2-D structure** (the publication-quality figure
Rfam uses). The pipeline is trim → sanitize → render:

```bash
bash scripts/run_r2r.sh
# Outputs land in work/15_msa_plots/<cluster>_r2r.pdf
#                   work/15_msa_plots/<cluster>_namemap.tsv
```

What the three internal steps do (transparent to the caller but worth
knowing for debugging):

- **Trim** (`trim_msa.py`): mlocarna's global alignment has 5′/3′ tails
  where only a minority of sequences have residues. The trimmer drops
  every column below `--min-coverage 0.5` (default), rebuilds SS_cons
  to keep only base pairs whose both ends survived.
- **Sanitize** (`sanitize_stockholm_for_r2r.py`): R2R 1.0.6 on Apple
  Silicon segfaults on `--GSC-weighted-consensus`, so we compute the
  required `#=GC cons` (majority IUPAC consensus with Rfam-style case),
  `#=GC conss` (visual-prominence tier), and `#=GC cov_SS_cons`
  (Watson-Crick pair frequency) annotations in pure Python and write a
  fully-annotated minimal Stockholm. Sequence names get renamed to
  `s001` / `s002` / … (R2R also crashes on names with `:/|.` or longer
  than ~25 chars); the original-name mapping lands in `<cluster>_namemap.tsv`.
- **Render**: R2R's normal render step takes the annotated Stockholm
  and emits PDF.

On a working R2R build (most Linux installs) you can collapse trim →
sanitize → render into the simpler trim → `r2r --GSC-weighted-consensus`
→ render flow; the hand-rolled annotations are only needed because of
the Apple Silicon segfault.

## 13. (Optional) Build CMs from the new alignments

Once the per-cluster alignments look right, you can build covariance
models with Infernal and use them to scan back across the library —
or any other library — for additional members:

```bash
mkdir -p work/16_cm
for stk in work/14_msa/*/results/result.stk; do
    cluster=$(basename "$(dirname "$(dirname "$stk")")")
    cmbuild --noss "work/16_cm/${cluster}.cm" "$stk"
    cmcalibrate "work/16_cm/${cluster}.cm"
done

# Re-search the pilot library
for cm in work/16_cm/*.cm; do
    name=$(basename "$cm" .cm)
    cmsearch --tblout "work/16_cm/${name}.tbl" "$cm" "$ref"
done
```

This closes the loop: any new hits from the new CMs can be fed back
into step 4 on the next iteration, which over time should grow the
catalogue of recognised xrRNA-like elements.

## Recap of where everything lands

| Stage | Output dir | Key artifact |
|---|---|---|
| Fragments | `work/05_fragments/` | per-BAM CSVs of 5′ ends |
| Calibration | `work/06_filter_v7/` | α̂ + Q-Q plot per replicate |
| Peak calls | `work/06_filter_v7/` | `{r}_peaks.csv`, `{r}_peaks.gtf` |
| CM hits | `work/08_infernal/` | `<cm>.tbl` |
| Combined peaks + seqs | `work/06_filter_v7/` | `123_peaks_combined_with_seqs.csv` (every row) |
| Annotated peaks | `work/10_postproc/` | annotated + per-accession CSVs (carries seq cols) |
| Downstream seqs (clustering FASTA) | `work/06_filter_v7/` | `*_downstream_seqs_unique.fa` (filtered) |
| Cluster-annotated peaks | `work/10_postproc/` | `123_peaks_annotated_with_clusters.csv` |
| IGV tracks | `work/11_tracks/` | bedgraphs + tile BEDs |
| Locarna scores | `work/06_filter_v7/seq_split/locarna_out/` | per-pair `.stk` + `.out` |
| Cluster FASTAs | `work/13_msa_groups_hybrid/` | one FASTA per cluster |
| MSAs | `work/14_msa/<cluster>/results/result.stk` | mlocarna output |
| MSA plots | `work/15_msa_plots/` | per-cluster PNG/SVG + R2R PDF |
| CMs | `work/16_cm/` | per-cluster `.cm` + `.tbl` |
