# viromelib + flaviLib — pipeline and figure reference — 2026-11-03

**Purpose.** A standalone index to the whole project: the environments and repositories, the analysis
pipeline, the most recent outputs and how to find them, and the figures with their inputs and
conventions. A person with access to the git repositories and to the laptop and cluster environments
should be able to locate the entire pipeline and its current outputs from this document alone.

**This document's location.** `~/Desktop/viromelib/workflows/FIGURES_log_consolidate_20261103.md`
(laptop, repo `viromelib-local`). `.md` documentation lives in `workflows/` by file kind (LAYOUT.txt 3).

**Authority.** The run logs are the record of what was run and are authoritative for commands, job IDs
and decisions. LAYOUT.txt is authoritative for where files go. Where this document and a source disagree,
section 13 lists it. Sources consolidated, with the date each is current to:

| document | home | current to |
|---|---|---|
| `README_261015-261019_re-run_viromelib.md` | cluster `$ROOT/workflows/` | 2026-10-21 |
| `README_261021_flaviLib_150x150.md` | cluster `$ROOT/workflows/` | 2026-10-28 |
| `FIGURES_status.md` | cluster `$ROOT/workflows/` | 2026-10-28 |
| `LAYOUT.txt` | cluster `$ROOT/LAYOUT.txt` | 2026-10-26 |
| `README.md` | cluster `$ROOT/README.md` | draft; section 13 |

---

## 1. Environments and repositories

### 1.1 The two machines

| | cluster (HPC, C2B2) | laptop |
|---|---|---|
| root | `$ROOT = /groups/as6282_gp/scratch_bkup/jgg2144/xrrna_seq/viromelib` | `~/Desktop/viromelib/` |
| role | reads, alignment, peak calling, library-wide `cmsearch`, annotation | CM building, structure work, R figures, IGV |
| conda env | `viromelib` | `viromelib` + R + R-scape 2.6.11 in `~/local/bin` |

Cluster login used in the run logs: `jgg2144@hpc.c2b2.columbia.edu`.

### 1.2 Three git repositories

| repo | remote | working tree | owns |
|---|---|---|---|
| `viromelib` (cluster) | `git@github.com:jgezelle/viromelib.git` | `$ROOT` | `scripts/`, `workflows/` run logs, `README.md`, `LAYOUT.txt`, CM models and alignments under `work/06_infernal/` (since 2026-10-26, commit `061683d`) |
| `viromelib-local` (laptop) | `git@github.com:jgezelle/viromelib-local.git` | `~/Desktop/viromelib/` | viromelib `.Rmd` notebooks, CM models under `work/infernal/cm/` |
| `flaviLib-local` (submodule) | `git@github.com:jgezelle/flaviLib-local.git` | `~/Desktop/viromelib/flaviLib/` | flaviLib notebooks and local analyses |

Submodule rule: commit inside `flaviLib/`, then commit the updated pointer in the parent (two commits).
Never `git add -A` in the parent (it would stage deletions of tracked CM models absent from the tree);
add explicit paths and read `git status --short` first. For scripts and run logs the cluster repo wins
over any differing copy; for notebooks the laptop copy is live while figures iterate.

### 1.3 Environment files — source one per shell, by absolute path

| file | run (`$DS`) / analysis (`$AN`) | reference / BED / taxonomy | role |
|---|---|---|---|
| `scripts/env_261015_re-run_viromelib.sh` | `261015_re-run_viromelib` / `fullBED` | 2,792-feature pre-tile / 8,178-tile / `virus_stats.csv` | frozen reproduction |
| `scripts/env_261019_orderedRef.sh` | `261019_re-run_viromelib_orderedRef` / `orderedRef` | 1,521-feature ordered / 5,967-tile / `virus_taxonomy_261019.csv` | **viromelib publication** |
| `scripts/env_261021_flaviLib.sh` | `261021_flaviLib_150x150` / `10mM`,`01mM` | flavi pre-tile / tile20 BED / `virus_taxonomy_flaviLib_261021.csv` | **flaviLib publication** |

`env_261019_orderedRef.sh` sets, among others:
`FA=pilot_library_20260319.uniqueIDs.ordered.fa`, `BED=pilot_library_20260319.uniqueIDs.ordered.bed`,
`FAMILY_CSV=virus_taxonomy_261019.csv`,
`LOG=$ROOT/work/00_logs/261019_re-run_viromelib_orderedRef`,
`POSTPROC_DIR=$ROOT/work/07_postproc/261019_re-run_viromelib_orderedRef/orderedRef`.
`$FAMILY_CSV` reaches step 11 through the environment; the job log prints `taxonomy: …` and silently
falls back to `virus_stats.csv` if the env was not sourced.

### 1.4 Raw data (cluster)

| what | where |
|---|---|
| original FASTQs (never modified) | `/groups/as6282_gp/data_bkup/jgg2144/<run folder>/` |
| viromelib pilot | `20260312_Aviti_ViromeLib_HighReads/` (Aviti #104, 150×0 SE) |
| flaviLib 150×150 | `20260903_Aviti_150x150_flaviLib-DMS/…/Samples/A####/` |
| links the pipeline reads | `$ROOT/work/01_raw/` (flat; symlinks) |
| sample sheet per run | `$ROOT/data/metadata/<run>_samples.tsv` (maps sample → X/R → replicate → link) |

### 1.5 File movement cluster → laptop

```
cluster  work/06_infernal/cm/*.cm        ◀─ built on laptop, scp up, md5 must match
cluster  work/06_infernal/<DS>/*.tbl     ◀─ step 10 cmsearch vs reference
cluster  work/07_postproc/<DS>/<AN>/     ─ scp + MANIFEST (md5) ▶  laptop figures/<fig>/inputs/
cluster  work/04_align, work/08_tracks   ─ scp ▶  laptop IGV inputs
```
Figure inputs come only from re-run outputs under `work/<stage>/<DS>/`, md5s recorded in each figure's
`inputs/MANIFEST.tsv`.

---

## 2. Libraries

### 2.1 viromelib

Aviti #104, 2026-03-12, single-end 150×0. 3 replicates × (X = Xrn1, R = RppH-only control); reps 1–2
prep r0369, rep 3 prep r0364. Tiles 196 nt at 80 nt step. Replicate gate ≥ 2 of 3.

| stage | accessions | viruses (TaxID) | features | oligos / tiles |
|---|---:|---:|---:|---:|
| GenBank records input to `parse.py` | 1,719 | — | — | — |
| designed pre-tile reference | 1,678 | 1,457 | 2,792 | 8,178 (visualization BED only) |
| **ordered library** (`library_v4.fa`) | **1,226** | **1,055** | **1,521** | **5,967** |
| ordered, Fig. 1 (Fiersviridae excluded) | 1,224 | 1,054 | — | — |

- Ordered pool: `data/refseqs/library_v4.fa`, md5 `c4ce9fc23ba31b0eb9c136b6195cb13c`.
- Ordered reference + BED: `data/refseqs/ordered_ref/pilot_library_20260319.uniqueIDs.ordered.{fa,bed}`,
  sequence-set md5 `a5b7d1c190908047fc3dd321d9dea522` (identical to RZ's independent build).
- Design → order gap: 1,268 features < 171 nt + 3 over `--max-ambig 4`. 0 unexplained.
- Fiersviridae (NC_055057.1, NC_074553.1; TaxID 2027243) are bacteriophages that passed the
  eukaryotic-host query through sample-level host annotation; reported in the NCBI-space figure, excluded
  from Fig. 1 (decision 23).

### 2.2 flaviLib (150×150 PE)

Aviti 2026-08-29, paired-end. Two Mg²⁺ conditions (10 mM, 1 mM), called separately, never pooled. 4
replicates per condition; reps 1–2 prep rJG0416, reps 3–4 prep rJG0441. Tiles 196 nt at 20 nt step.
Replicate gate ≥ 2 of 4 (decision 26 open). 439 accessions, 273 TaxIDs (ratio 1.61). 28 accessions with
no family rank; `AB698917.1` is a dsRNA virus (inclusion unconfirmed).

---

## 3. Pipeline

Driver (cluster): `bash $ROOT/scripts/run_viromelib_pipeline.sh [FROM_STEP]` runs steps 5–12 as one
SLURM `afterok` chain. Each wrapper reads the sample sheet and writes an `inputs.md5` beside its outputs.

| step | wrapper / script | output under `work/<stage>/<DS>/[<AN>/]` |
|---|---|---|
| 5 align | `align_bowtie2_se.sbatch <RUN> <REF>` | `04_align/` BAM + index |
| 6 fragments | `make_fragments_se.sbatch <RUN>` | `05_fragments/` |
| 7 calibrate | `calibrate_null.sbatch <RUN> <AN> <TILE_BED>` | `06_filter/<AN>/calibration/alphas.tsv` |
| 8 peaks | `call_peaks_v73.sbatch <RUN> <AN> <TILE_BED> <MIN_SF> <MAX_Q>` | `06_filter/<AN>/rep<r>_peaks.{csv,gtf}` |
| 9 combine | `combine_reps.sbatch <RUN> <AN> <FASTA>` | `06_filter/<AN>/peaks_combined_with_seqs.csv` |
| 10 CM search | `cmsearch_library.sbatch <RUN> <REF> <CM_DIR> <CM…>` | `06_infernal/<DS>/<CM>.cm.tbl` (per reference) |
| 11 annotate | `postprocess_peaks.sbatch <RUN> <AN> <NAME> <CM…>` | `07_postproc/<AN>/<NAME>.csv` |
| 12 tracks | `make_tracks.sbatch <RUN> <REF>` | `08_tracks/` bedGraphs |
| 13 final sites | `apply_publication_cutoffs.py` (+ `stall_distance_to_stop.py`, `flag_shared_oligos.py`, `make_release_tracks.py sites`) | `07_postproc/<AN>/*.sites_final.csv`, BED, GTF, xlsx |
| 17 coverage | `featureCounts` → `coverage_tables.py` | `07_postproc/<AN>/coverage.*.tsv` |
| 18 NCBI space | `ncbi_virus_space.py`; 18b `family_absence_reasons.py` | figure inputs |
| 19 IGV | `make_igv_figures.py` | IGV sessions + PNG/SVG |

v7.3 peak caller is used unchanged; the only pipeline input that differs between runs is the tile BED.
Scripts live in `$ROOT/scripts/` (cluster repo). Full per-script notes: `README.md` and LAYOUT.txt 4.

---

## 4. Runs and their most recent outputs

| library | run (`$DS`) / analysis | status | feeds |
|---|---|---|---|
| viromelib | `261015_re-run_viromelib` / `fullBED` | frozen reproduction of the 2026-07-21 table | no figure |
| viromelib | `261019_re-run_viromelib_orderedRef` / `orderedRef` | **publication** | every viromelib figure |
| flaviLib | `261021_flaviLib_150x150` / `10mM`,`01mM` | **publication** | every flaviLib figure |
| flaviLib | 260807 75×75 | retired (`vintage_75x75/`) | no figure |

### 4.1 viromelib orderedRef — latest outputs on the cluster (verified 2026-11-03)

Postproc directory `$POSTPROC_DIR = $ROOT/work/07_postproc/261019_re-run_viromelib_orderedRef/orderedRef/`:

| file | md5 | notes |
|---|---|---|
| `peaks_annotated.csv` | `4f452542f85c9c185a10bc75f258ea5a` | 2026-10-19 20:43, 24,187 rows, **5 CMs** (1a, 1b, 2, cl3, **cl4** — the pre-v9 Class 4). No dumbbell |
| `peaks_annotated_per_accession.csv` | `f285cc4967455d28b9f424f03f1d6244` | collapsed view |

`peaks_annotated_inputs.md5` records the inputs that annotation used:

| input | md5 |
|---|---|
| `work/06_filter/261019_re-run_viromelib_orderedRef/orderedRef/peaks_combined_with_seqs.csv` | `1892587819f8b73aac1b4550412c5c1b` |
| `data/virus_taxonomy_261019.csv` | `a6f566b42fcb193dfa1ece145e665137` |
| `work/06_infernal/261019_re-run_viromelib_orderedRef/CL_Class_1a_xrRNA_final.cm.tbl` | `5770619a4b790eb2ef0653fc192737b8` |
| `.../CL_Class_1b_xrRNA_final.cm.tbl` | `5dd488936b671a6ce3380f060276829c` |
| `.../CL_Class_2_xrRNA_final.cm.tbl` | `ca6faddd777b3856c49cb29eed1aba4f` |
| `.../cl3_NAR_2025.…SS_cons.cm.tbl` | `377de72dedadd56f1edd1984b0c94f87` |
| `.../cl4.cm.tbl` | `337c560ebf8866bce745537a8d74efc7` |
| `scripts/post_process_peaks.py` | `42e6d86c9aafe5dd008dd6ee48e347b8` |

CM search tables: `$ROOT/work/06_infernal/261019_re-run_viromelib_orderedRef/` (five `.tbl`, listed above).
Final sites: `$POSTPROC_DIR/261019_re-run_viromelib_orderedRef_orderedRef.sites_final.csv`.

**Current annotation carries five CMs and the pre-v9 Class 4; it has no dumbbell and no class4_v9.**
Figure 1 requires the six known models, so a re-annotation is required — section 11.

Result numbers (unchanged by re-annotation, which does not touch the gates):

| | value |
|---|---|
| candidate rows | 24,187 |
| tile rows passing ≥ 2 of 3 | 446 |
| **unique sites** | **391** (233 in all 3 replicates) |
| features / accessions / viruses with a site | 301 / 286 / 264 |
| sites with a CM hit (five-CM annotation) | 83 |

Text phrasing: "391 stall sites on 264 viruses (286 accessions), all v7.3 gates in ≥ 2 of 3 replicates;
233 in all three."

### 4.2 flaviLib 261021 — latest outputs

Per-condition postproc `$ROOT/work/07_postproc/261021_flaviLib_150x150/<COND>/`:
`261021_flaviLib_150x150_<COND>.peaks_annotated.csv` (seven CMs), sites and `summary.json` from step 13a.
Superseded annotation kept as the record: `work/07_postproc/260903_trueRef_final/`.

| | 10 mM | 1 mM |
|---|---:|---:|
| candidate rows | 34,102 | 40,324 |
| unique sites (≥ 2 of 4) | 699 | 459 |
| in all 4 replicates | 585 | 386 |
| accessions with a site | 315 | 245 |
| sites with a CM hit (7 CMs) | 330 | 243 |

Detectability floor (decision 29): 123 control reads, best tile; marginal to 2,312; bands 9 / 28 / 402 of
439. Mg²⁺ pairing: 784 sites either condition (337 10 mM only, 111 1 mM only, 336 both); median paired
ΔSF −0.029 vs unpaired 0.588 / 0.404 — the effect is loss of detection.

### 4.3 NCBI virus space (viromelib supplementary)

1,007 of 20,969 +ssRNA eukaryotic-host species (4.80 %); 53 of 68 families.

---

## 5. CM models

### 5.1 The set

`cm_set_261021` (decision 3, 2026-10-21), seven models:
`CL_Class_1a_xrRNA_final`, `CL_Class_1b_xrRNA_final`, `CL_Class_2_xrRNA_final`,
`cl3_NAR_2025.align.mpkadded.edited.confident.dedup.IDcurated.SS_cons`, `class4_cl4_v9_20260824`,
`class1c`, `dumbbell_RF00525`.

| location | holds |
|---|---|
| laptop `work/infernal/cm/cm_set_261021/` | the 7 models (CM building happens in `work/infernal/`) |
| cluster `work/06_infernal/cm/` (`$CMDIR`) | the flat applied set; superseded models in `cm/vintage/` (including the pre-v9 `cl4.cm`) |

All 7 md5s verified end to end across laptop, cluster `$CMDIR`, and the run-log checksums on 2026-10-26.
Re-verify after any model changes, before re-annotating. The hand-off (laptop → `$CMDIR`) is manual; the
md5 is the guard.

### 5.2 Known before this study vs developed in it

Decided 2026-11-03. Figure 1, and any "with a CM match" count, use only models known going in:

| model | status | shown in Fig. 1 |
|---|---|---|
| Class 1a, 1b, 2 (`CL_Class_*_xrRNA_final`) | pre-existing xrRNA classes | yes |
| Class 3 (`cl3_NAR_2025…`) | pre-existing (NAR 2025) | yes |
| Class 4 (`class4_cl4_v9_20260824`) | class pre-existing (Dilweg 2023); **v9 model refit in-study** (decision 9) | yes, v9 model |
| Dumbbell (`dumbbell_RF00525`) | published Rfam RF00525; a 3′UTR element, not an xrRNA class — kept as a possible resistant hit | yes, labelled "Dumbbell" |
| `class1c` | defined in this study (pestivirus 1c; decision 33 lineage) | no |
| expanded class 1b, orthoflavivirus subclasses | to be defined in-study | no |

So Figure 1's annotation uses **six CMs**: 1a, 1b, 2, cl3, class4_v9, dumbbell. `class1c` and future
in-study models stay out of the Figure-1 annotation. Legends that mix shown models with in-study ones
elsewhere must state which were known beforehand (decision 34). The panel title may change from "CM
matches" to "resistant hits", since dumbbell is included as a resistant element rather than an xrRNA.

---

## 6. Figure layout (laptop)

Per-figure directory:
```
<lib>/figures/<figure>/
  <notebook>.Rmd            live notebook
  archive/                  superseded notebooks
  inputs/                   FLAT; run or library in the FILENAME; MANIFEST.tsv of md5s
  plots/<YYYYMMDD>_<strategy>/   outputs; never overwritten
```
Input filenames carry the run, e.g. `261019_re-run_viromelib_orderedRef.peaks_annotated_known6cm.csv`,
`261021_flaviLib_150x150_01mM.peaks_annotated.csv`.

| directory | notebook (live) | library |
|---|---|---|
| `figures/fig1/` | `fig1A_v8_4panel.Rmd` (v7 → `archive/`) | viromelib |
| `figures/rank/` | `rank_plot_v9_20261019.Rmd` | viromelib |
| `figures/ncbi_virus_space/` | `ncbi_virus_space_20261020.Rmd` | viromelib |
| `figures/representation/` | `representation_floor_20261026.Rmd` (one notebook, `LIB` switch, two figures) | both |
| `figures/IGV/` | from `scripts/make_igv_figures.py` (cluster) | both |
| `flaviLib/figures/rank/` | `rank_plot_flaviLib_v1_20261021.Rmd` | flaviLib |
| `flaviLib/figures/input_distribution/` | `inputdist_flaviLib_v1_20261021.Rmd` | flaviLib |
| `flaviLib/figures/mg_dependence/` | `mg_dependence_flaviLib_20261026.Rmd` | flaviLib |
| (unplaced) | `sf_vs_depth_20261027.Rmd` | both |
| `figures/coverage/` | not started | both |
| `figures/dms/`, `figures/sbcnv_5end/` | self-contained sub-projects, own `scripts/`+`workflows/` | viromelib |

CM class distribution is not a separate figure: two panels inside the fig1 notebook (outputs; ordered and
detected inputs).

---

## 7. Inputs per figure and their cluster source

Each input is copied from a re-run output, renamed with its run prefix, and recorded in
`inputs/MANIFEST.tsv`.

### 7.1 viromelib Fig. 1A (`figures/fig1/inputs/`)

| laptop name | cluster source |
|---|---|
| `261019_re-run_viromelib_orderedRef.peaks_annotated_known6cm.csv` | the six-CM re-annotation (section 11) under `$POSTPROC_DIR` |
| `virus_taxonomy_261019.csv` | `$ROOT/data/virus_taxonomy_261019.csv` (md5 `a6f566b4…`) |
| `pilot_library_20260319.uniqueIDs.ordered.fa` | `$ROOT/data/refseqs/ordered_ref/` |
| `accession_presence_by_sample_261019.csv` | cluster figure inputs |
| `ncbi_family_counts_261019.csv` | cluster figure inputs (Entrez 2026-10-19) |

### 7.2 viromelib rank plot (`figures/rank/inputs/`)

`261019_re-run_viromelib_orderedRef.peaks_annotated_known6cm.csv` (same file as 7.1);
`261019_re-run_viromelib_orderedRef_orderedRef.sites_final.csv` (`$POSTPROC_DIR`);
`virus_taxonomy_261019.csv`.

### 7.3 flaviLib figures

`261021_flaviLib_150x150_<COND>.peaks_annotated.csv`, the step-13a sites files,
`virus_taxonomy_flaviLib_261021.csv`, and `flaviLib_input_representation.tsv` /
`flaviLib_detectability_floor.txt` from `detectability_floor.py --prefix flaviLib_`.

---

## 8. Plotting conventions

| topic | convention |
|---|---|
| greys | `GREY_REFSEQ` (`grey88`) = in RefSeq, not in the library (panel 1 only). `GREY_LIB` (`grey45`) = in the library; the denominator of every panel, including "XR-positive, no hit". No other use of either |
| models shown | only models known before the study (section 5.2): classes 1a, 1b, 2, 3, 4 and dumbbell. `class1c` and future in-study models are not drawn and do not count toward "with a CM match" |
| dumbbell | labelled "Dumbbell" (not a numbered class); last in single-class priority; colour `#8C564B`; a resistant element, not an xrRNA |
| excluded taxa | Fiersviridae excluded from viromelib Fig. 1 by family and accession, with a `stopifnot` guard; stated in the caption |
| names | `organism` (GenBank ORGANISM) for labels; `species` (ICTV binomial) only as the formal name. NC_034223.1 = "Jingmen tick virus (MGTV/V4/11)" |
| units | state accessions or viruses (TaxIDs) with every count |
| denominator | the ordered library, never the designed reference |
| XR (stall) definition | treated SF ≥ 0.3, control SF ≤ 0.1, fold change ≥ 3, stall reads ≥ 10, q ≤ 0.01, in ≥ MIN_REPS replicates |
| `peaks_annotated.csv` | holds every candidate row regardless of gates; reapply `stall_reads ≥ 10` when plotting an "all" set |
| IGV | shared locked axis, CPM; axis label says CPM |
| supplementary numbering | never in file or directory names |
| no-family group | shown as a labelled group or explicitly dropped |

---

## 9. Figure status

| figure | status | needs before final |
|---|---|---|
| viromelib Fig. 1A | v8 built 2026-11-03, not knit | six-CM re-annotation (section 11); decision 32 (denominator) |
| viromelib rank plot | needs change | path update (11.4); HIGHLIGHT name; `stall_reads ≥ 10` on `all`; bands after decision 32; six-CM columns |
| flaviLib rank plot | needs change | remove HCV; negative control (decision 30); `stall_reads ≥ 10`; bands; resolve `rank_matches_input_stall_site_rank: false` |
| flaviLib input distribution | fixes unverified | a knit |
| NCBI virus space | done | — |
| representation + floor | flaviLib ready (floor 123) | viromelib floor (decision 32) |
| Mg²⁺ dependence | built, unverified | knit; lead with the detection-loss panel |
| SF and p vs depth | partly knit (1 mM) | 10 mM; place in a directory |
| CM class distribution (in fig1) | not started | decisions 28, 35, 36; input-vs-output structure of `cm_distributions_viromeLib_flaviviridae.Rmd` |
| IGV | viromelib done | flaviLib not started |
| coverage | not started | — |

### 9.1 Fig. 1A v8 (2026-11-03)

1. Greys per section 8, with legend keys for both.
2. Fiersviridae excluded at every input, by family and accession, with a `stopifnot` guard.
3. Caption panel order matches drawn order (NCBI, XR-positive, CM, lengths).
4. Only models known before the study are drawn (section 5.2): classes 1a, 1b, 2, 3, 4 and dumbbell;
   `class1c` and future in-study models excluded; `any_hit` recomputed over the shown models, with the
   file's value kept as `any_hit_all_models`.
5. Class 4 read from `class4_cl4_v9_20260824.cm_hit` when present (`CLASS4_PREFERRED`), else `cl4.cm_hit`.
6. `PEAKS_FILE` points at the six-CM re-annotation; the `stopifnot` guard halts the knit if a shown
   model's column is absent, so the figure cannot ship missing a known class. The knit prints the CM
   columns found and records the peaks md5 and CM columns in the QC log.
7. Output to `figures/fig1/plots/20261103_v8_greyfix_noFiers/`.

---

## 10. Open decisions affecting figures

| # | decision | affects |
|---|---|---|
| 2 | viromelib replicate cutoff ≥ 2 of 3 (current tables use it) | every viromelib count |
| 4 | site merge distance 0 vs 3 nt (viromelib 385 → 377 on fullBED; flaviLib 66 pairs within 1 nt) | every site count |
| 26 | flaviLib gate ≥ 2 of 4 vs ≥ 3 of 4 | flaviLib counts |
| 27 | flaviLib rank plot per condition or combined | flaviLib rank plot |
| 28 | which CM columns a figure shows; `class1c`/dumbbell mapping | CM panels, rank plots (resolved for Fig. 1 by section 5.2) |
| 30 | flaviLib negative control | flaviLib rank plot |
| 31 | tile-count/depth gradient: length or sub-pool | methods, representation figure |
| 32 | viromelib detectability floor | Fig. 1 denominator, rank plot bands |
| 35 | `cmstat` on CMs (class4_v9 first) | CM panels |
| 36 | minus-strand cmsearch hits (`--toponly`) | CM panels |
| 37 | PLLG "CM-blind" re-test on full features | claims a clade lacks a known class |
| — | flaviLib `rank_matches_input_stall_site_rank: false` | flaviLib rank plot ordering |

Settled and applied: 3, 9, 18, 19, 20, 21, 22, 23, 24, 25, 29, 33, 34, and (2026-11-03) the known-going-in
rule for Fig. 1.

---

## 11. Commands — six-CM re-annotation of viromelib orderedRef

The gates do not use CM annotation, so the site set stays at 391; only the CM columns, `any_hit` and the
CM-hit count change. Confirm wrapper argument order against the run log before submitting. New outputs
carry a `known6cm` suffix; the five-CM 2026-10-19 files stay as the record.

### 11.1 Confirm the two new models exist on the cluster and match the laptop set

```bash
source /groups/as6282_gp/scratch_bkup/jgg2144/xrrna_seq/viromelib/scripts/env_261019_orderedRef.sh
CMDIR=$ROOT/work/06_infernal/cm
md5sum $CMDIR/class4_cl4_v9_20260824.cm $CMDIR/dumbbell_RF00525.cm | tee $LOG/known6cm_cm_md5.txt
# these two md5s must equal the laptop cm_set_261021 and the run-log checksums
```

### 11.2 cmsearch the two models missing for the ordered reference

The ordered-reference `.tbl` set holds only the five pre-v9 tables (section 4.1). `class4_cl4_v9` and
`dumbbell_RF00525` must be searched against the reference; 1a, 1b, 2, cl3 are reused unchanged.

```bash
INF=$ROOT/work/06_infernal/$DS; mkdir -p $INF/logs
for CM in class4_cl4_v9_20260824 dumbbell_RF00525; do
  [ -f $INF/$CM.cm.tbl ] && { echo "exists: $CM"; continue; }
  sbatch -o $INF/logs/cmsearch_%j.out $ROOT/scripts/cmsearch_library.sbatch \
    $DS $ROOT/data/refseqs/ordered_ref/$FA $CMDIR $CM.cm
done
# when done: strand check (decision 36; tblout col 10)
for t in $INF/class4_cl4_v9_20260824.cm.tbl $INF/dumbbell_RF00525.cm.tbl; do
  printf '%-40s plus %4s minus %4s\n' "$(basename $t)" \
    "$(grep -v '^#' $t | awk '$10=="+"' | wc -l)" "$(grep -v '^#' $t | awk '$10=="-"' | wc -l)"
done | tee $LOG/known6cm_strand.txt
```

### 11.3 Annotate with the six known models (step 11)

```bash
sbatch -o $POSTPROC_DIR/logs/postproc_known6cm_%j.out $ROOT/scripts/postprocess_peaks.sbatch \
  $DS $AN peaks_annotated_known6cm \
  CL_Class_1a_xrRNA_final.cm CL_Class_1b_xrRNA_final.cm CL_Class_2_xrRNA_final.cm \
  cl3_NAR_2025.align.mpkadded.edited.confident.dedup.IDcurated.SS_cons.cm \
  class4_cl4_v9_20260824.cm dumbbell_RF00525.cm
# verify
cat $POSTPROC_DIR/logs/postproc_known6cm_*.out          # must print taxonomy: .../virus_taxonomy_261019.csv
F=$POSTPROC_DIR/peaks_annotated_known6cm.csv
md5sum $F; head -1 $F | tr ',' '\n' | grep '\.cm_hit$'   # expect 6 columns, incl class4_cl4_v9 and dumbbell
cat $POSTPROC_DIR/peaks_annotated_known6cm_inputs.md5
```

### 11.4 Final sites (step 13), then confirm the site set is unchanged

Re-run the 2026-10-19 step-13 calls with the new input and a `known6cm` prefix. Recover the exact
command lines from the shell history:

```bash
history | grep -E 'apply_publication_cutoffs|stall_distance_to_stop|flag_shared_oligos|make_release_tracks' | tail -20
# compare old vs new site set and CM-hit count
python3 - $POSTPROC_DIR/261019_re-run_viromelib_orderedRef_orderedRef.sites_final.csv \
          $POSTPROC_DIR/261019_re-run_viromelib_orderedRef_orderedRef.known6cm.sites_final.csv <<'PY'
import sys, pandas as pd
a,b=(pd.read_csv(f) for f in sys.argv[1:3]); k=["accession","stall_site_genomic_position"]
print("sites old/new:",len(a),len(b),"| same set:",set(map(tuple,a[k].values))==set(map(tuple,b[k].values)))
print("cm_hit cols new:",[c for c in b.columns if c.endswith(".cm_hit")])
PY
```

Expect 391 sites in both and an identical site set.

### 11.5 Download to the laptop and point the notebook at it

```bash
SRC=/groups/as6282_gp/scratch_bkup/jgg2144/xrrna_seq/viromelib/work/07_postproc/261019_re-run_viromelib_orderedRef/orderedRef/peaks_annotated_known6cm.csv
F1=~/Desktop/viromelib/figures/fig1
scp jgg2144@hpc.c2b2.columbia.edu:"$SRC" \
    $F1/inputs/261019_re-run_viromelib_orderedRef.peaks_annotated_known6cm.csv
( cd $F1/inputs && md5 -r *.csv *.fa > MANIFEST.tsv && cat MANIFEST.tsv )
```

`fig1A_v8_4panel.Rmd` already points `PEAKS_FILE` at
`261019_re-run_viromelib_orderedRef.peaks_annotated_known6cm.csv` and reads Class 4 from the v9 column.
No notebook edit is needed once the file is in place.

### 11.6 Directory flattening and the rank plot

`fig1/inputs/` still has the legacy `20261019_rerun_viromelib_ordered/` subdirectory. Flatten it (move the
four non-peaks inputs up, drop the superseded five-CM `peaks_annotated.csv`), give `figures/rank/` its own
`inputs/`, and update the rank notebook's paths and HIGHLIGHT name (NC_034223.1 =
"Jingmen tick virus (MGTV/V4/11)"). Commit with explicit paths:

```bash
cd ~/Desktop/viromelib
git add workflows/FIGURES_log_consolidate_20261103.md \
        figures/fig1/fig1A_v8_4panel.Rmd figures/fig1/archive/fig1A_v7_4panel.Rmd \
        figures/fig1/inputs/MANIFEST.tsv \
        figures/rank/rank_plot_v9_20261019.Rmd figures/rank/inputs/MANIFEST.tsv
git status --short
git commit -m "Figures: consolidated reference 20261103; fig1A v8 (6 known CMs incl dumbbell)"
```

---

## 12. Quick index — where to find things

| looking for | path |
|---|---|
| what was run, job IDs, decisions | cluster `$ROOT/workflows/README_*.md` |
| where any file goes | cluster `$ROOT/LAYOUT.txt` |
| viromelib publication peaks | `$ROOT/work/07_postproc/261019_re-run_viromelib_orderedRef/orderedRef/` |
| viromelib CM search tables | `$ROOT/work/06_infernal/261019_re-run_viromelib_orderedRef/` |
| flaviLib publication peaks | `$ROOT/work/07_postproc/261021_flaviLib_150x150/<COND>/` |
| CM models (applied) | cluster `$ROOT/work/06_infernal/cm/` ; (built) laptop `work/infernal/cm/cm_set_261021/` |
| references and BEDs | `$ROOT/data/refseqs/ordered_ref/` |
| taxonomy tables | `$ROOT/data/virus_taxonomy_261019.csv`, `virus_taxonomy_flaviLib_261021.csv` |
| figure notebooks | laptop `~/Desktop/viromelib/<lib>/figures/<figure>/` (section 6) |
| scripts | cluster `$ROOT/scripts/` |

---

## 13. Statements in the source documents that are out of date

None changes a result.

| where | says | current |
|---|---|---|
| `README.md` | run log `README_261015_re-run_viromelib.md` | `README_261015-261019_re-run_viromelib.md` |
| `README.md` | `inputs/<run>/MANIFEST.tsv` | flat `inputs/`, run in filename |
| `README.md` | `scripts/virus_stats_ranked.py` | `scripts/virus_taxonomy.py` |
| `README.md` | submodule repo `jgezelle/flaviLib` | `jgezelle/flaviLib-local` |
| `README.md` | ordered pool `library_v4_1.nospace.fa` | `library_v4.fa` |
| `FIGURES_status.md` | floor "must be re-derived on best-tile" | done: 123 (decision 29) |
| `FIGURES_status.md` | gradient Spearman +0.688 | +0.745 (best-tile) |
| `FIGURES_status.md` | notebooks under `local_READMEs/` | live copies under each `figures/<figure>/` |
| viromelib log 12e | flaviLib floor 122 median | 123 best-tile |
| viromelib log, progress table | steps 13–15 pending | 13 done both runs; 14, 15 not run |
| viromelib log 19c | rank plot HIGHLIGHT already renamed | the copy on hand still says "Mogiana tick virus" |
| flaviLib log, decision 35 | two statuses merged in one cell | class1c done; class4 next; five pending |
| `fig1A_v7` caption | lengths = panel 2 | drawn fourth (fixed in v8) |
