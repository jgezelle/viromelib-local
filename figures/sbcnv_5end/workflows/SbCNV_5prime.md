# SbCNV-5 5′-End Mapping: Running Log

**Goal:** Map the 5′ ends of a pool of RNAs that were isolated and sequenced together, against the SbCNV-5 reference.
**Started:** October 9, 2026
**Project root:** `/groups/as6282_gp/scratch_bkup/jgg2144/xrrna_seq/viromelib`
**Reference:** `$ROOT/data/refseqs/nc_024077_sbcnv-5.fasta` (NC_024077)

```bash
ROOT=/groups/as6282_gp/scratch_bkup/jgg2144/xrrna_seq/viromelib
REF="$ROOT/data/refseqs/nc_024077_sbcnv-5.fasta"
```

## Summary (Oct 9, 2026)

- **Sample:** A3480, 5′-ligation library (Illumina adapter + N8), gene-specific RT primer annealing at 18780–18798; 163,800 read pairs; no Xrn1/RppH/PNK treatment, so only native 5′-monophosphate ends were captured.
- **Pipeline:** trim N8 + adapters (and DMS-library crosstalk) → QC → align R1 end-to-end to a strain-patched NC_024077 (44 substitutions, same coordinates) → count 5′ ends per nt → figure (R) and IGV tracks.
- **Main 5′ ends** (130,098 on-target reads): **18681** (29.0%; element 2 +1), **18142** (8.2%; the ~800 bp gel band, 657-nt product), 18682 (7.6%), **18471** (3.5%), **18576** (3.0%; element 1 +2), 18680 (3.0%), **18279** (1.1%), 18574/18575 (element 1 start, +1), 18743 (0.9%).
- **Elements 1 and 2** show the same pattern: ends at element start, +1 and +2, plus a cluster 3–7 nt upstream.
- **Assay stall sites** 18574, 18679 and 18680 coincide with sequenced 5′ ends; the strongest end is 1–2 nt downstream of each. Stall sites past 18798 cannot be measured with this RT primer.
- **Longest products:** 18116 (683 nt, ≥ 10 reads), 18047 (752 nt, ≥ 3 reads). All major ends lie in the 3′ UTR.
- **Background:** ~33% of reads, a smooth rise toward the primer (length bias) with minor real ends on top.

## 1. Sample

| ID | Link name | Well | i7 | i5 | Read pairs | R1+R2 |
|---|---|---|---|---|---|---|
| A3480 | 260930-scnv-blue-5end-A3480 | E1 | oSL1287 TGTGACTG | oSL1231 ACGTTCAG | 163,800 | 29.4 MB |

- Run: Element Aviti 150x150, `20260929_AV100007_092826_Sternberg_Aviti`, Bases2FastQ `093026_Bases2FastQ`.
- Depth is ~22% of the 750k target; judged usable at import.
- Download, QC and linking are documented in `README_20260930.md` (job `835671`, Oct 7).

## 2. Inputs

```
$ROOT/work/01_raw/
├── 260930-scnv-blue-5end-A3480_R1.fastq.gz   -> symlink into data_bkup run dir
└── 260930-scnv-blue-5end-A3480_R2.fastq.gz
```

Raw data source: `/groups/as6282_gp/data_bkup/jgg2144/20260930_Aviti_150x150_DMS-virome-scnv/ngs_data/.../Samples/A3480/`

**Reference features (from J.G., Oct 9):** genome 19,199 nt. Positions 18574–18760 hold two near-repeats; by functional-element boundaries:

| Element | Coordinates (1-based) | Sequence |
|---|---|---|
| Element 1 (upstream) | 18574–18644 | `cccgtcagcttgacgttaactgccactttggttgaagtgtgatcaaccgtgcctggggcgagcatatcggc` |
| Element 2 (downstream) | 18680–18750 | `cccgtcagcctgacgttaactgccactttggttgaagtgtgatcagccatgcccggagcgagcatatcggt` |

The two elements are 106 nt apart (start to start) and 71 nt each. Coordinates confirmed by exact match against `$REF`.

**Genome annotation (J.G., Oct 9):** EW45_gp1 CDS (polyprotein) **351–18068**; inferred 3′ UTR **18069–19199**. All major 5′ ends lie in the 3′ UTR (18142 is 74 nt into it); only a few reads have ends in the last ~150 nt of the CDS.

## 3. Session setup

Paths live in `$ROOT/scripts/config_sbcnv5p.sh` (exports `ROOT`, `REF`, `RAW`, `TRIM`, `QC`, `ALN`, `IDX`, `THREADS`, `N8`; creates output dirs).

```bash
# login node
ROOT=/groups/as6282_gp/scratch_bkup/jgg2144/xrrna_seq/viromelib
srun --pty -t 04:00:00 -c 8 --mem=16G --chdir="$ROOT" /bin/bash -l

# inside the job
conda activate viromelib
source scripts/config_sbcnv5p.sh
for t in cutadapt bowtie2 samtools bedtools fastqc multiqc; do printf "%-10s " $t; command -v $t || echo MISSING; done
```

## 4. Layout

`work/` already holds output from earlier viromelib analyses (`02_trim` … `07_postproc`, `06_filter_v*`).
To keep this analysis separate, each step writes into a `sbcnv5p/` subfolder:

```
work/
├── 01_raw/                       # shared links (inputs)
├── 02_trim/sbcnv5p/              # N8- and adapter-trimmed reads
├── 03_qc/sbcnv5p/                # FastQC / MultiQC, trimming reports
└── 04_align/sbcnv5p/             # BAMs, 5′-end counts, tracks
    └── patched/                  # step 04b: patched reference + final counts (used downstream)
```

Code and metadata in the repo:

```
$ROOT/
├── scripts/
│   ├── config_sbcnv5p.sh          # paths (source first)
│   ├── 02_trim_sbcnv5p.sh
│   ├── 03_qc_sbcnv5p.sh
│   ├── 04_align_5p_sbcnv5p.sh
│   ├── 04b_patch_ref_sbcnv5p.sh
│   └── 05_plot_5p_sbcnv5p.R        # figure (run locally, §6.1)
├── workflows/
│   ├── SbCNV_5prime.md             # this README
│   └── fig_sbcnv_5end.Rmd          # local R Markdown driver for the figure (§6.1)
└── data/metadata/
    └── sbcnv5_annotations.tsv      # CDS, UTR, elements, RT primer, assay stall sites
```

## 5. Library structure

5′ adapter ligation. The ligated adapter is the Illumina sequencing adapter, with **8 random nt (N8) at its 3′ end**, so N8 sits directly against the RNA 5′ end:

```
5′-[Illumina adapter]-NNNNNNNN-[RNA 5′ end ...]-3′
                      └─ N8 ─┘ └─ first templated base = 5′ end to call
```

Expected 5′-end read: **R1**, beginning `N8` then RNA base 1 (to confirm by peeking at reads, §7 Q2).

**Handling N8:** hard-trim exactly 8 nt from the start of R1 (cutadapt `-u 8`) instead of letting the aligner soft-clip it:
- Soft-clipping is unreliable at an exact boundary. By chance the last N matches the reference base just upstream of the true 5′ end ~25% of the time (the last two ~6%), so the aligner extends into it and calls the end 1 to 2 nt too far upstream. That smears every real end into a false stutter.
- After a fixed 8-nt trim, alignment is end-to-end, so R1 start = RNA 5′ end exactly.
- **Not used as a UMI / no deduplication** (decided Oct 9; not needed for this analysis). Counts are raw read 5′ ends.

**Adapters** (same TruSeq-style i5/i7 as the DMS libraries, `xrRNASHAPE/DMS` step 02):
- R1 3′ read-through: i7 `AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC`
- R2 3′ read-through: i5 rc `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`
- When R2 reads through, the 8 nt before the adapter are rc(N8) and remain on R2's 3′ end. Removing them conditionally isn't clean in cutadapt, and R2 isn't needed to call 5′ ends, so 5′ ends are called from **R1 alone**. R1 has no random bases on its 3′ side.

**RT (gene-specific primer, Oct 9):** RT primer `CTGGAGTTCAGACGTGTGCTCTTCCGATCT` + `cgtctggtgacgagctaat`
- The 3′ part anneals to `attagctcgtcaccagacg` = **+ strand 18780–18798** (unique; no other site within 4 mismatches).
- The 5′ tail is the i7-side adapter, which is why R1 reads straight from the insert into `AGATCGGAAGAGCACACGTCTG…` (§9.3).
- So every on-target cDNA runs from the RNA 5′ end to **18798**: product length = 18798 − 5′ end + 1. The 18798 "3′ end" is the primer site, not a biological end. R2 always starts with the primer; R2 reads that don't are off-target priming (e.g. rRNA).
- Bases 18780–18798 in reads come from the primer, so any variant there isn't informative. (The fixed 18779 T→C is 1 nt upstream of the primer, so it is real RNA sequence.)

**Length bias:** RT processivity, PCR and clustering all favor short products, so peak heights are not proportional to RNA abundance across very different product lengths (e.g. element 2 at ~118 nt vs element 1 at ~225 nt). Within a peak cluster, relative heights are more comparable.

**Multimapping:** elements 1 and 2 differ at 6 / 71 positions. Because on-target R1 reads always run from the 5′ end into unique downstream sequence (up to 142 nt), the closest alternative placement for any R1 ≥ 18 nt is ≥ 10 mismatches away. Best-hit alignment is safe; MAPQ is still reported.

## 6. Plan

| Step | What | Tool | Status |
|---|---|---|---|
| 01 peek | Confirm N8 is at the start of R1 (base composition, raw reads) | zcat, awk | done Oct 9 (§9.1) |
| 02 trim | `scripts/02_trim_sbcnv5p.sh`: pass 0 drop DMS `LEAD5` crosstalk; pass 1 `-u 8` on R1, i7/i5 read-through, `-q 20`, `-m 18` | cutadapt | done Oct 9 (§9.5) |
| 03 qc | `scripts/03_qc_sbcnv5p.sh`: FastQC raw + trimmed, MultiQC (incl. cutadapt logs); pair survival through pass 0/1; trimmed R1 pos 1–4 composition, R1 lengths, % R2 on-target | FastQC, MultiQC, awk | done Oct 9 (§9.6) |
| 04 align + ends | `scripts/04_align_5p_sbcnv5p.sh`: R1 single-end, best hit, `--end-to-end --very-sensitive`; tag on-target reads (R2 = RT primer); 5′ end per read; counts per position (MAPQ ≥ 10); bedGraphs; summary incl. longest product | bowtie2, samtools, bedtools | done Oct 9 vs NC_024077 (§9.6) |
| 04b patch | `scripts/04b_patch_ref_sbcnv5p.sh`: 2-round substitution-only consensus of NC_024077 from the R1 alignments (coordinates unchanged), then re-run step 04 into `04_align/sbcnv5p/patched/` | samtools consensus (≥1.16), bowtie2 | done Oct 9 (§9.7); **patched counts are the ones used downstream** |
| 04 (cont.) | Summary also bins on-target product lengths (50 nt) with estimated library size (product + 144 bp) to compare with the gel | awk | added Oct 9 |
| 05 figure | `scripts/05_plot_5p_sbcnv5p.R` + `data/metadata/sbcnv5_annotations.tsv`: A genome map, B 5′ ends 17,900–19,199 (log reads), C elements 1–2 at 1-nt resolution; tracks for CDS, elements, assay stalls, RT primer, strain variants. Writes .svg, .pdf, `_linear_break.svg/.pdf`, `_peaks.tsv`, `_stalls.tsv` | R: ggplot2, patchwork, ggrepel, svglite | v2 Oct 9: CDS/UTR + assay stall sites; unmeasured zone shaded. Run locally (§6.1) |

**Order:** 02 trim → 03 qc → 04 align + ends → 05 figure.

```bash
bash scripts/02_trim_sbcnv5p.sh
bash scripts/03_qc_sbcnv5p.sh
bash scripts/04_align_5p_sbcnv5p.sh
bash scripts/04b_patch_ref_sbcnv5p.sh
```

**Key design points**
- Which read starts at the RNA 5′ end determines everything downstream. If R1 begins at the RNA 5′ end, the 5′ end is the R1 start (leftmost base on + strand, rightmost on − strand). If the library is antisense to the RNA, R2 carries it instead.
- End-to-end alignment is required: soft-clipping the read start would shift called 5′ ends. Non-templated bases at the 5′ end (e.g. from the ligation or RT step) should be reported, not hidden.
- Low depth (163,800 pairs) is plenty for a single small genome, but positions with a handful of reads should be treated with caution.

### 6.1 Local analysis: figure and IGV (laptop)

Steps 02–04b run on the cluster. The figure and IGV browsing run locally from a few small files.

**Local layout** (`~/Desktop/viromelib/figures/fig_sbcnv_5end/`):

```
fig_sbcnv_5end/
├── 05_plot_5p_sbcnv5p.R                 # same script as $ROOT/scripts/
├── workflows/fig_sbcnv_5end.Rmd         # same file as $ROOT/workflows/
└── 261009/                              # one folder per run (date)
    ├── inputs/
    │   ├── 260930-scnv-blue-5end-A3480.5p_counts.tsv   # from work/04_align/sbcnv5p/patched/
    │   ├── substitutions.tsv                           # from work/04_align/sbcnv5p/patched/
    │   └── sbcnv5_annotations.tsv                      # from data/metadata/
    └── plots/                           # written by the Rmd
        ├── 261009_SbCNV5_5prime_ends.svg / .pdf
        ├── 261009_SbCNV5_5prime_ends_linear_break.svg / .pdf  # linear axis with a scale break
        ├── 261009_SbCNV5_5prime_ends_peaks.tsv         # every 5′ end: reads, %, product length, nearest feature
        └── 261009_SbCNV5_5prime_ends_stalls.tsv        # each stall site: reads at/±1/±2, strongest end within ±5 nt, offset
```

Copy the inputs from the cluster (replace `<cluster>` with the login address):

```bash
D=~/Desktop/viromelib/figures/fig_sbcnv_5end/261009/inputs; mkdir -p $D
P=/groups/as6282_gp/scratch_bkup/jgg2144/xrrna_seq/viromelib
scp jgg2144@<cluster>:$P/work/04_align/sbcnv5p/patched/{260930-scnv-blue-5end-A3480.5p_counts.tsv,substitutions.tsv} $D/
scp jgg2144@<cluster>:$P/data/metadata/sbcnv5_annotations.tsv $D/
```

**Running the figure:**
1. R packages: ggplot2 (≥ 3.4), patchwork, ggrepel, svglite, knitr. To install only what's missing:
   ```r
   pkgs <- c("ggplot2", "patchwork", "ggrepel", "svglite", "knitr")
   missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
   if (length(missing)) install.packages(missing)
   ```
2. Open `workflows/fig_sbcnv_5end.Rmd` in RStudio and Knit. The setup chunk sets the working directory for every chunk with `knitr::opts_knit$set(root.dir = ...)` (a plain `setwd()` only lasts one chunk). The settings chunk points the script at `261009/inputs/` and `261009/plots/` via `options(sbcnv5.*)`; for a new run, change `run <- "261009"`.
3. Command-line alternative (from `fig_sbcnv_5end/`; `--counts` also accepts the 5′-end bedGraph):
   ```bash
   Rscript 05_plot_5p_sbcnv5p.R --counts 261009/inputs/260930-scnv-blue-5end-A3480.5p_counts.tsv \
     --subs 261009/inputs/substitutions.tsv --ann 261009/inputs/sbcnv5_annotations.tsv \
     --out 261009/plots/261009_SbCNV5_5prime_ends.svg
   ```
4. Edit the SVG in Illustrator/Inkscape for final layout. Windows, label thresholds and colors are in the settings block at the top of the R script; features and stall sites are rows in `sbcnv5_annotations.tsv`.

**From reads to the plot (how each stem is computed):**

1. **One 5′ end per read** (step 04/04b, `…5p_reads.tsv.gz`). Each trimmed R1 starts at the RNA 5′ end, because the N8 was hard-trimmed and alignment is end-to-end (no soft-clipping). For a read on the + strand the 5′ end is its leftmost aligned base (BAM `POS`); on the − strand it would be the rightmost base. The row also records MAPQ and whether the read is **on-target** (its R2 starts with the RT primer `CGTCTGGTGACGAGCT`).
2. **Counts per nucleotide** (`…5p_counts.tsv`): reads with MAPQ ≥ 10 are tallied by (position, strand). `n_all` = all such reads; `n_ontarget` = on-target only.
3. **5′-end bedGraph** (`…5p_ontarget.plus.bedGraph`, and `A3480.5p_ontarget.counts.bedGraph` in the IGV bundle): the `n_ontarget` column for the + strand, one line per nucleotide with ≥ 1 read. bedGraph is 0-based half-open, so `18680 18681 37780` means 37,780 5′ ends at nt 18681. It is the same data as step 2, in a browser format. The `.pct` version divides by the 130,098 on-target reads with ends ≤ 18760.
4. **The figure** reads either file (`--counts` accepts the counts table or the bedGraph; both give identical results, checked Oct 9). It keeps + strand, on-target, `pos ≤ 18760` (ends closer to the primer are mostly primer sequence), and draws **one stem per nucleotide whose height is the number of reads with their 5′ end there**, on a log axis starting just below 1 so single-read positions stay visible. Labeled peaks: in B, positions ≥ 0.8% of on-target reads, one label per 10-nt cluster; in C, every position ≥ 1%. Product length in the tables = 18798 − 5′ end + 1.

Nothing is smoothed, binned or normalized in the stems: each one is a raw count at a single nucleotide.

**What the figure shows:** on-target reads only (R2 starts with the RT primer), + strand, patched-reference counts at MAPQ ≥ 10, 5′ ends ≤ 18760 (ends closer to the primer are mostly primer sequence). A: genome map (CDS, 3′ UTR, zoom box). B: 5′ ends per nt, 17,900–19,199, log read axis; region past 18760 shaded as not measured; tracks for CDS/UTR, elements, assay stall sites (filled = high score), RT primer, strain variants. C: elements 1–2 at single-nt resolution.

**Second figure, `_linear_break.svg`:** the same counts on a **linear** read axis, so the size differences between the dominant ends read directly (the log axis compresses them). The y axis is broken once: a lower segment from 0 to just above the second-largest end (0–12,000) and an upper segment around the largest (35,000–42,000, holding only 18681). Both segments use the same reads-per-mm scale; panel heights are proportional to their ranges, and slashes on the y axis mark the break. A: 17,900–18,760 with ends ≥ 1,000 reads labeled. B, C: element 1 and element 2 starts as single-nt bars, read counts above bars and positions on the axis. Ranges are set automatically from the data (`lo_max`, `hi_rng` in the script).

**IGV** (`sbcnv5_igv_bundle.zip`, unzip anywhere, e.g. `figures/IGV/`):
- `sbcnv5_patched.fa` (+ `.fai`): strain sequence in NC_024077 coordinates, loaded as the genome.
- `A3480.5p_ontarget.counts.bedGraph`: 5′-end reads per nt; `A3480.5p_ontarget.pct.bedGraph`: % of on-target reads. bedGraph is 0-based: `18680 18681 37780` = 37,780 ends at nt 18681.
- `sbcnv5_features.bed` (CDS, UTR, elements, primer, stall sites), `sbcnv5_strain_variants.bed`.
- File → Open Session → `sbcnv5_5prime_session.xml` (bar charts, no windowing, counts on log scale, locus 17,900–18,800). File → Save SVG Image exports the current view.

## 7. Open questions

1. ~~Library structure~~ Answered Oct 9: 5′ adapter ligation, Illumina adapter with 3′ N8 (§5).
2. ~~Which read = RNA 5′ end~~ Confirmed Oct 9: **R1**, N8 at positions 1–8 (§9).
3. ~~3′ side of the insert~~ Answered Oct 9: gene-specific RT primer annealing at 18780–18798 (§5).
4. ~~RNA 5′ end chemistry~~ Answered Oct 9: no Xrn1, RppH or PNK treatment. Only native 5′-monophosphate ends were ligatable, so every captured end (major peaks and background) was already 5′-P in the RNA; the major peaks are interpreted as stable 5′-P degradation products.
5. **Expected ends:** any known or predicted 5′ ends to check against?
6. ~~R2 fixed start~~ Resolved Oct 9: it is the 3′ end of the dominant 118-nt insert itself (genomic, ~18798), read directly from the i7 side with no extra handle (§9.3). Whether that 3′ end is the RNA's own end or an RT-primer site is Q3.
7. ~~Adapter family~~ Answered Oct 9: TruSeq-style i5/i7, same as DMS libraries (§5).
8. ~~DMS crosstalk~~ Confirmed Oct 9: `LEAD5` `GGATTAATATAAT` is not in the reference (either strand), so R1 reads starting with it are crosstalk from the DMS libraries; pass 0 removes them.
9. ~~Reference mismatches / repeats~~ Resolved Oct 9: the repeats are elements 1 and 2 (§2). The dominant insert maps uniquely to element 2 (5 mismatches vs 28 at element 1) and carries 5 fixed differences from NC_024077 plus 2 polymorphic sites (§9.4). Fixed differences likely strain-specific (J.G., Oct 9); align to NC_024077, best hit.
10. ~~~800 bp gel band~~ Resolved Oct 9: real product. 5′ end **18142** (10,603 reads, 11.3% of on-target, product 657 nt → ~801 bp library), as predicted (~18,140). The 651–700 nt bin holds 12.2% of on-target reads. Not a heteroduplex artifact.
11. ~~Low MAPQ from strain divergence~~ Resolved Oct 9 by step 04b: 44 substitutions; MAPQ < 10 21.1% → 2.1%, aligned 83.1% → 93.0% (§9.7). Losses had been uneven (short products hit hardest).
12. **Primer-only products:** R1 reads of products < ~40 nt lie mostly inside the 19-nt primer site (18780–18798), which is primer-derived, not RNA. Exclude 5′ ends > 18760 (product < 39 nt, < 20 nt of templated RNA) from peak calls and the figure.
13. ~~Annotation for the figure~~ Answered Oct 9: EW45_gp1 CDS (polyprotein) 351–18068; inferred 3′ UTR 18069–19199; assay stall sites (§9.8).

## 8. Log

- **Oct 9, 2026:** README started. Inputs confirmed in `work/01_raw` (R1/R2 links for A3480). Plan drafted.
- **Oct 9, 2026:** Library structure recorded (5′ ligation, Illumina adapter + N8). Decided to hard-trim N8 (fixed 8 nt) rather than rely on soft-clipping. R1 expected to carry the 5′ end. N8 not used as a UMI; no dedup. Next: peek at reads to confirm N8 on R1; find out what's on the 3′ side.
- **Oct 9, 2026:** Added `scripts/config_sbcnv5p.sh` and session setup (§3): interactive `srun` job, `conda activate viromelib`, tool check. Next: run the read peek (§6 step 01) inside the job.
- **Oct 9, 2026 11:29:** Interactive job on `c1205` (8 cores, 16G, 4 h). `viromelib` env has cutadapt, bowtie2, samtools, bedtools, fastqc, multiqc (all in `~/.conda/envs/viromelib/bin`). `bedGraphToBigWig` not yet checked (only needed for tracks).
- **Oct 9, 2026 11:30:** Read peek (first 100k reads). R1 confirmed as 5′-end read with N8. R2 starts with a fixed sequence in ~90% of reads. Details in §9.
- **Oct 9, 2026 11:43:** Adapters taken from the DMS trimming step (TruSeq-style i5/i7). Spotted DMS `LEAD5` handle at the start of some R1 reads (crosstalk). Wrote `scripts/02_trim_sbcnv5p.sh` (pass 0 crosstalk filter, pass 1 N8 + adapters). Decided to call 5′ ends from R1 alone.
- **Oct 9, 2026 11:47:** Reference check (exact 25-mer / 13-mer, both strands). Dominant R1 species starts at **+18681**. R2 fixed sequence: no hit. `LEAD5`: no hit, so crosstalk confirmed. Trim on hold until R2 sequence is resolved.
- **Oct 9, 2026 11:49:** R2 fixed sequence resolved: it is the 3′ end of the dominant insert (118 nt, ~18681–18798), read straight into the i7 adapter with no handle. 52.8% of R1 reads contain it. Trim script needs no change. Possible tandem repeat / reference mismatches in this region (§7 Q9).
- **Oct 9, 2026 11:56:** Got the reference FASTA and element boundaries (element 1 18574–18644, element 2 18680–18750). Dominant 5′ end 18681 is 1 nt inside element 2. The insert differs from NC_024077 at 5 fixed + 2 polymorphic positions (§9.4).
- **Oct 9, 2026 12:02:** RT was gene-specific (anneals 18780–18798), so all products end at 18798 and product length = 18798 − 5′ end + 1. Decided: best-hit alignment to NC_024077 (multimapping between elements not expected, checked). Question to answer: is the element 2 peak dominant because of RT/PCR length bias or abundance; what is the longest product. Wrote `scripts/04_align_5p_sbcnv5p.sh`. Goal: figure of all 5′ peaks on the genome (step 05).
- **Oct 9, 2026 12:17:** J.G. saw a major ~800 bp library band on the gel; predicted 5′ end ~18,140 (adapter-corrected). Added product-length bins + estimated library size to step 04 summary. Order set: trim → qc → align → figure.
- **Oct 9, 2026 12:20:** Gel detail: faint ~200 bp band (likely the 118-nt product, expected ~262 bp), darker higher band, 600–900 bp not resolved. TapeStation/Bioanalyzer would size it properly.
- **Oct 9, 2026 12:25:** Wrote `scripts/03_qc_sbcnv5p.sh` (modeled on the DMS QC block).
- **Oct 9, 2026 12:30:** Step 02 run. Crosstalk 3,594 pairs (2.2%) removed; 151,259 pairs kept (92.3% of raw). 57,189 R1 reads untrimmed at 142 nt = inserts ≥ 142 nt (long products or off-target); to be placed by alignment.
- **Oct 9, 2026 12:38:** Steps 03 and 04 run (vs NC_024077). QC clean (N8 gone; 90.1% on-target). The ~800 bp band is real: 5′ end 18142 (657-nt product). Element 1 shows the same end pattern as element 2. 21% of aligned reads lost at MAPQ < 10 to strain divergence → wrote step 04b (patched reference) and re-run before the figure.
- **Oct 9, 2026 12:46:** Step 04b run: 44 substitutions vs NC_024077 (no deletions). Aligned 93.0%, MAPQ < 10 2.1%, mean NM 1.03. Peak positions unchanged; short products recovered most (51–100 nt bin 1.4% → 9.8%). Patched counts adopted. Next: figure (step 05).
- **Oct 9, 2026 13:45:** Wrote the R figure script and annotation table; first draft rendered from the patched counts (130,098 on-target reads at 663 positions, ends ≤ 18760). Waiting on CDS and assay stall-site coordinates.
- **Oct 9, 2026 13:58:** Added CDS/UTR and 7 assay stall sites to the figure; stall comparison in §9.8. Sites at 18574/18679/18680 coincide with sequenced 5′ ends (dominant end 1–2 nt downstream); 18831/18928/18929/19089 lie downstream of the RT primer and cannot be measured in this library.
- **Oct 9, 2026 14:02:** Checked user coordinates against the sequence (all consistent, §9.9). Characterized background 5′ ends (32.6% of reads, rising toward the primer). Built IGV bundle (strain FASTA, feature/stall BED, variant BED, session XML).
- **Oct 9, 2026 14:16:** No Xrn1/RppH/PNK treatment: all captured ends are native 5′-P. Figure route settled: R script for the publication figure, IGV for exploration.
- **Oct 9, 2026 14:38:** Added `fig_sbcnv_5end.Rmd` for running the figure locally (`~/Desktop/viromelib/figures/fig_sbcnv_5end/`, inputs in `261009/`, outputs in `plots/`). R script now reads settings from `options(sbcnv5.*)` and creates the output folder. Test-knitted OK.
- **Oct 9, 2026 15:04:** Added summary and §6.1 (local figure + IGV workflow). Rmd/R script to be copied to `$ROOT/workflows/` and `$ROOT/scripts/`; annotations to `data/metadata/`.
- **Oct 9, 2026 15:22:** README §6.1 now explains how each stem is computed from reads → counts → bedGraph. Plot script accepts the 5′-end bedGraph as input (identical output to the counts table).
- **Oct 9, 2026 15:28:** Added a second figure (`_linear_break`): linear read axis with one scale break (0–12k / 35k–42k), overview + element 1/2 zooms.
- **Oct 9, 2026 15:38:** Repo is under git; this README lives at `workflows/SbCNV_5prime.md`. Commit plan: download/link scripts, SbCNV-5 5′ pipeline, and unrelated changes as separate commits.

## 9. Results

### 9.1 Read peek (Oct 9, first 100k reads)

**R1:** positions 1–8 are close to uniform (each base 0.15–0.34), then composition becomes strongly skewed from position 9 (C 0.50 at pos 9). Confirms N8 at the start of R1, and R1 position 9 = RNA 5′ end.
- Mild G depletion in the N8 (0.15–0.18), a common synthesis/ligation bias; not a concern.
- One species dominates: 6 of the first 10 reads are N8 + `CCGTCAGCCTGACGTTAACTGCCACTTTGGTTGAAGTGCAATCAGCC...`.
- Its 5′ end looks heterogeneous by a few nt: most start `CCGTCAG`, one `CGTCAG` (1 nt shorter), one `CCCGTCAG` (with 3 extra nt before the N8 boundary). This is exactly the case the fixed 8-nt trim handles; the aligned 5′ ends will show the real distribution.
- Other reads include rRNA (e.g. `...GTACACACCGCCCGTCACACCATGGGAG...`, the conserved SSU rRNA 3′ region), expected as background.

**R2:** no random bases. ~90% of reads begin with one fixed sequence, `CGTCTGGTGACGAGCTAATGATCACGCCGTCGCTACTCGCC...` (composition ≥0.89 for one base at every position 1–16). The rest are other templates (incl. rRNA). Origin open (§7 Q6).

### 9.2 Reference check (Oct 9, exact match on linearized `$REF`)

| Query | Source | + strand (1-based start) | − strand |
|---|---|---|---|
| `CCGTCAGCCTGACGTTAACTGCCAC` | dominant R1 species (after N8) | **18681** | none |
| `CGTCTGGTGACGAGCTAATGATCAC` | fixed start of ~90% of R2 | none | none |
| `GGATTAATATAAT` | DMS `LEAD5` handle | none | none |

- The dominant R1 species is sense (+ strand) with its main 5′ end at **18681** (reads beginning `CCGTCAG`); the `CGTCAG` and `CCCGTCAG` variants in §9.1 would be 18682 and ~18680. Alignment will give the real distribution.
- The R2 fixed sequence (≥60 nt identical across reads in §9.1) is not in the reference as an exact 25-mer.

### 9.3 R2 fixed sequence and insert structure (Oct 9)

**12-mer tiling of the R2 consensus** (`CGTCTGGTGACGAGCTAATGATCACGCCGTCGCTACTCGCCGTGATGGACCGATGTGCTC`):

| R2 pos | 12-mer | + | − |
|---|---|---|---|
| 1 | CGTCTGGTGACG | | 18787 |
| 7 | GTGACGAGCTAA | | 18781 |
| 13, 19, 25, 31 | | | none |
| 37 | TCGCCGTGATGG | | 18539, 18645, 18751 |
| 43 | TGATGGACCGAT | | 18533, 18745 |
| 49 | ACCGATGTGCTC | | 18527 |

- R2 position 1 is the rc of + strand **18798**: the R2 fixed start is the 3′ end of the sense insert, not a handle.
- Positions 13–36 have no exact 12-mer hit: the RNA differs from the reference there (or the reference is off).
- Positions 37–60 hit sites spaced ~106 nt apart: looks like a tandem repeat in the reference.

**R1 read-through:** 86,557 / 163,800 R1 reads (52.8%) contain rc(R2 start) `CATTAGCTCGTCACCAGACG`. Example (N8 removed):

```
CCGTCAGCCTGACGTTAACTGCCACTTTGGTTGAAGTGCAATCAGCCATGCCCGAAGCGAGCACATCGGTTCATCACGGCGAGTAGCGACGGCGTGATCATTAGCTCGTCACCAGACG | AGATCGGAAGAGCACACGTCTGAA
└──────────────────────────────── 118-nt insert, 5′ end 18681 ───────────────────────────────────────────────┘   └ i7 adapter
```

The insert runs straight into the i7 adapter, so the trim script needs no extra handle. Reads also vary at two internal positions (`CCCGAAG`/`CCCGGAG`, `GGTTCATC`/`GGTCCATC`).

**Preview of 5′ ends** (assumes the 3′ end is fixed at 18798; 5′ end = 18788 − R1 position of the match):

| R1 pos of match | Reads | Insert (nt) | Implied 5′ end |
|---|---|---|---|
| 107 | 39,152 | 118 | **18681** |
| 106 | 10,343 | 117 | 18682 |
| 108 | 4,081 | 119 | 18680 |
| 113 | 2,397 | 124 | 18675 |
| 111 | 1,667 | 122 | 18677 |
| 112 | 1,129 | 123 | 18676 |
| 45 | 1,197 | 56 | 18743 |
| 1 | 1,228 | (no N8) | artifact |

Alignment will replace this estimate.

### 9.4 Dominant insert vs reference (Oct 9)

Ungapped comparison of the 118-nt insert (two read variants from §9.3) against `$REF`:
- Best hit **18681–18798** (5 mismatches); next best is the element 1 copy at 18575 (28 mismatches). Unique placement in element 2.
- The 5′ end 18681 is **1 nt inside element 2** (element 2 starts at 18680; the minor 18680 end is exactly the element boundary). The minor 18675–18677 cluster sits just upstream of element 2.

| Position | Ref | Reads | Where |
|---|---|---|---|
| 18719 | T | C | element 2 (fixed in both variants) |
| 18720 | G | A | element 2 (fixed) |
| 18735 | G | A or G | element 2 (polymorphic) |
| 18744 | T | C | element 2 (fixed) |
| 18751 | C | T or C | just past element 2 (polymorphic) |
| 18764 | C | T | downstream (fixed) |
| 18779 | T | C | downstream (fixed) |

All substitutions, no indels, so NC_024077 coordinates still hold. For alignment: 5–7 mismatches fit within bowtie2's default end-to-end score threshold at 118 nt, but short reads (e.g. the 56-nt 18743 product) carry 3–4 of them and are near the limit. Option: a patched reference with the 5 fixed substitutions (same coordinates).

### 9.5 Trimming (step 02, Oct 9)

| Stage | Pairs | % raw |
|---|---|---|
| Raw | 163,800 | 100 |
| After pass 0 (`LEAD5` crosstalk removed: 3,594) | 160,206 | 97.8 |
| After pass 1 (N8 + adapters; 8,947 too short < 18 nt) | 151,259 | 92.3 |

Adapter found: R1 100,873 (63.0%), R2 95,722 (59.7%) of pass-1 input.

**Trimmed R1 lengths** (top 15). For on-target reads, length < 142 = full product, and 5′ end = 18798 − length + 1:

| R1 length | Reads | Implied 5′ end | Note |
|---|---|---|---|
| 142 | 57,189 | ≤ 18657 | untrimmed: product ≥ 142 nt (long products, incl. element 1 and the ~800 bp band, or off-target) |
| 118 | 39,812 | 18681 | dominant |
| 117 | 10,821 | 18682 | |
| 119 | 4,220 | 18680 | element 2 boundary |
| 124 / 122 / 123 / 121 / 126 / 125 / 120 | 2,468 / 1,711 / 1,176 / 947 / 916 / 795 / 754 | 18675 / 18677 / 18676 / 18678 / 18673 / 18674 / 18679 | upstream cluster |
| 141 | 1,341 | (18658) | likely 142-nt reads with 1 nt trimmed (partial adapter / quality) |
| 56 | 1,317 | 18743 | |
| 75 | 794 | 18724 | |
| 80 | 766 | 18719 | |

Matches the §9.3 preview. Over a third of kept R1 reads (57,189) come from products ≥ 142 nt; alignment will show where their 5′ ends lie.

### 9.6 QC and alignment vs NC_024077 (steps 03, 04; Oct 9)

**QC:** trimmed R1 position 1 is C 0.52 (skewed, N8 fully removed). R2 starts with the RT primer in 136,306 / 151,259 pairs (**90.1% on-target**). MultiQC: `work/03_qc/sbcnv5p/260930-scnv-blue-5end-A3480_multiqc.html`.

**bowtie2 (R1, end-to-end, NC_024077):** 151,259 reads; 83.06% aligned (79.37% once, 3.69% >1 time); 25,616 unaligned. Of 125,643 aligned: 118,910 on-target (94.6%); **26,544 (21.1%) MAPQ < 10**. NM mode 5–6 (strain differences), tail to ≥ 12.
- Since only 3.69% multimap, most MAPQ < 10 reads are unique alignments with many mismatches (bowtie2 lowers MAPQ when the score is near its threshold). With ~14,950 off-target pairs, roughly 10k on-target reads also failed to align. Both losses are likely uneven across the genome → step 04b.

**Top on-target 5′ ends (+ strand, MAPQ ≥ 10; % of on-target counted):**

| 5′ end | Reads | % | Product (nt) | Context |
|---|---|---|---|---|
| 18681 | 29,512 | 31.4 | 118 | element 2 start +1 |
| **18142** | **10,603** | **11.3** | **657** | ~800 bp gel band |
| 18682 | 7,475 | 8.0 | 117 | element 2 +2 |
| 18471 | 4,464 | 4.8 | 328 | |
| 18576 | 3,914 | 4.2 | 223 | element 1 start +2 |
| 18680 | 3,323 | 3.5 | 119 | element 2 start |
| 18675 / 18677 / 18676 / 18673 / 18678 / 18674 / 18679 | 1,967 / 1,335 / 824 / 619 / 588 / 540 / 375 | 6.5 total | 120–126 | 1–7 nt upstream of element 2 |
| 18279 | 1,469 | 1.6 | 520 | |
| 18575 / 18574 | 1,384 / 1,335 | 1.5 / 1.4 | 224 / 225 | element 1 start +1 / start |
| 18569 / 18571 / 18570 | 895 / 823 / 537 | 2.4 total | 228–230 | 3–5 nt upstream of element 1 |
| 18646 / 18645 / 18647 | 761 / 599 / 339 | 1.8 total | 152–154 | just past element 1 end (18644) |
| 18619 / 18453 / 18599 | 499 / 495 / 438 | | 180 / 346 / 200 | |

- **Elements 1 and 2 show the same pattern:** ends at element start, +1, +2, plus a cluster a few nt upstream (element 2: +1 dominant; element 1: +2 dominant).
- **Longest products:** most upstream end with ≥ 3 reads = **18116** (683 nt); a single read at 15379 (3,420 nt) is likely spurious.

**Product-length bins (on-target, + strand):** 101–150 nt 53.6% (~269 bp lib); 201–250 nt 12.8%; **651–700 nt 12.2% (~819 bp lib)**; 151–200 nt 7.0%; 301–350 nt 6.4%; 501–550 nt 2.4%; 401–450 nt 2.0%; everything else < 1.5% each.

### 9.7 Patched reference (step 04b, Oct 9) — counts used from here on

Consensus: round 1 38 substitutions, round 2 **44** (no deletions); NC_024077 coordinates unchanged. List: `04_align/sbcnv5p/patched/substitutions.tsv`.

| | NC_024077 | Patched |
|---|---|---|
| Aligned R1 | 125,643 (83.06%) | **140,599 (92.95%)** |
| Unaligned | 25,616 | 10,660 |
| Aligned >1 time | 3.69% | 5.48% |
| MAPQ < 10 | 26,544 (21.13%) | **2,949 (2.10%)** |
| Mean NM | 5.59 | **1.03** |
| On-target | 118,910 | 132,968 |

**Top on-target 5′ ends (patched, MAPQ ≥ 10):**

| 5′ end | Reads | % | Product (nt) | vs NC_024077 |
|---|---|---|---|---|
| 18681 | 37,780 | 29.0 | 118 | +28% |
| 18142 | 10,629 | 8.2 | 657 | unchanged |
| 18682 | 9,925 | 7.6 | 117 | +33% |
| 18471 | 4,541 | 3.5 | 328 | +2% |
| 18576 | 3,950 | 3.0 | 223 | +1% |
| 18680 | 3,890 | 3.0 | 119 | +17% |
| 18675 / 18677 / 18676 / 18678 / 18673 / 18674 / 18679 | 2,301 / 1,679 / 1,120 / 905 / 797 / 717 / 686 | 6.2 total | 120–126 | |
| 18279 | 1,470 | 1.1 | 520 | unchanged |
| 18575 / 18574 | 1,411 / 1,359 | 1.1 / 1.0 | 224 / 225 | |
| **18743** | 1,169 | 0.9 | 56 | newly above MAPQ 10 |
| 18569 / 18571 / 18570 | 930 / 860 / 568 | 1.8 total | 228–230 | |
| 18646 / 18645 | 795 / 634 | 1.1 total | 153–154 | |
| **18724** / **18738** | 710 / 578 | 0.5 / 0.4 | 75 / 61 | newly above MAPQ 10 |
| 18453 | 583 | 0.4 | 346 | |

- Peak positions are the same as before; only heights changed. Short products gained most (they carry the dense strain differences at 18719–18779 relative to their length), long products barely changed. So the NC_024077 run under-counted short products.
- Product-length bins: 101–150 nt 49.7%; 51–100 nt **9.8%** (was 1.4%); 201–250 nt 9.8%; 651–700 nt 8.9%; 151–200 nt 6.4%; 301–350 nt 5.5%; 251–300 nt 2.8%; 401–450 nt 2.2%; 501–550 nt 1.7%; 1–50 nt 1.3% (see §7 Q12).
- Longest products: ≥ 10 reads **18116** (683 nt); ≥ 3 reads **18047** (752 nt); single read at 15379 ignored. Products > 700 nt: 62 reads total.

### 9.8 Assay stall sites vs sequenced 5′ ends (Oct 9)

Stall sites from J.G.'s degradation assay (high score = 18574, 18680, 18928). From `SbCNV5_5prime_ends_stalls.tsv` (patched counts, on-target):

| Stall | High score | Measurable here | Reads at site | ±1 nt | ±2 nt | Strongest end within ±5 nt | Offset |
|---|---|---|---|---|---|---|---|
| 18574 | yes | yes | 1,359 | 2,989 | 7,270 | 18576 (3,950) | +2 |
| 18679 | | yes | 686 | 5,481 | 44,940 | 18681 (37,780) | +2 |
| 18680 | yes | yes | 3,890 | 42,356 | 53,186 | 18681 (37,780) | +1 |
| 18831 | | no (downstream of RT primer) | | | | | |
| 18928 | yes | no | | | | | |
| 18929 | | no | | | | | |
| 19089 | | no | | | | | |

- All three measurable stall sites are real 5′ ends in the sequencing, at both element starts. The strongest sequenced end sits **1–2 nt downstream** of each assay site, consistently. Worth checking whether that reflects the assay's position convention (e.g. last protected nt vs first nt of the product) or a real 1–2 nt trimming.
- The four sites past 18798 cannot be tested with this library (the gene-specific RT primer copies only RNA upstream of 18798). A primer closer to the genome 3′ end would be needed.
- Major sequenced ends with no assay stall nearby: **18142, 18279, 18471**, and 18743.

### 9.9 Coordinate checks, background, IGV (Oct 9)

**Coordinates checked against the sequence:**
- CDS 351–18068: starts ATG at 351; 17,718 nt (in frame); single in-frame stop at 18066–18068 (`TAA` in NC_024077, `TGA` in this strain via the 18067 A→G variant, still a stop).
- Elements 1/2, RT primer site: exact matches (§2, §5).
- 5′ end convention: reads beginning `CCGTCAG…` map to 18681, and the reference at 18681 reads `CCGTCAG…`, so a 5′ end coordinate is the first templated nucleotide.
- Stall sites: 18574 and 18680 are the same position in the two repeats, the first C of `…ttag|CCCGTCAG…`; 18679 is the G just upstream. The dominant sequenced ends are the 2nd and 3rd C (18575/18576, 18681/18682).
- Observation: the 18142 and 18279 ends are both followed by `TGGAGG` (`TATGGAGG…`, `GTTGGAGG…`).

**Background (everything outside ±2 nt of the 12 main peaks):** 42,367 reads (32.6%) over 624 positions, median 22 reads per position. The 11 positions ≥ 1% hold 60.7% of reads; the 123 positions ≥ 0.1% hold 88.2%.
- Density rises smoothly toward the primer: 0.3 reads/nt in the CDS end, 9.5 at 18069–18200, 22 at 18201–18400, 61 at 18401–18573, 107–175 across elements 1–2.
- Weak sequence preference: background 5′ nt is G-depleted (18% vs 31% in the region), slightly enriched for A/C/T.
- Not RT drop-off: the 5′ adapter is ligated to the RNA before RT, so a cDNA that stops early never copies the adapter and cannot be amplified. Background reads are ligated RNA 5′ ends: a broad population of decay/processing intermediates, inflated toward the primer by the length bias of RT, PCR and clustering. The smooth gradient is the length bias; the jagged structure on top (e.g. 18187–18201, 18379–18390, 18505–18513, 18537–18544, 18564–18572, 18618–18622, 18637, 18645–18650) is minor real ends.
- Ways to pin it down: a library from an intact in vitro transcript of the 3′ UTR (handling background + length bias), and knowing the 5′ end chemistry (§7 Q4).

**IGV bundle** (`sbcnv5_igv_bundle.zip`): `sbcnv5_patched.fa` (+ .fai; strain sequence, NC_024077 coordinates), `sbcnv5_features.bed` (CDS, UTR, elements, primer, stall sites), `sbcnv5_strain_variants.bed`, `sbcnv5_5prime_session.xml` (bar chart, log scale, no windowing, locus 17,900–18,800). Now also contains the 5′-end bedGraphs, so no copy from the cluster is needed:
- `A3480.5p_ontarget.counts.bedGraph`: reads per nt whose 5′ end is at that position (same values as `patched/…5p_ontarget.plus.bedGraph`, with an IGV track line).
- `A3480.5p_ontarget.pct.bedGraph`: frequency, % of the 130,098 on-target reads with 5′ ends ≤ 18760 (18681 = 29.04%).
bedGraph coordinates are 0-based half-open: the line `18680  18681  37780` is the 5′ end at nt 18681.
