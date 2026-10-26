# IGV sessions for top hits 

#### Aug 27, 2026
local download of necessary input files
### Received 11_tracks contents from RZ, copied to work dir as `11_tracks`
```bash
ROOT=~/Desktop/viromelib/
TRACKS=~/Desktop/viromelib/work/11_tracks
SCRIPTS=~Desktop/viromelib/scripts

```
```bash
conda activate viromelib 
cd $ROOT

    python scripts/make_igv_session.py \\
        --genome data/refseqs/pilot_library_20260319.uniqueIDs.fa \\
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
        --bed work/data/refseqs/pilot_library_20260319.uniqueIDs.bed \\
        --out 3replicates.xml
```
-------

#### Aug 25, 2026
```bash 
conda activate viromelib

ROOT=~/Desktop/viromelib
IGVDIR="$ROOT/figures/IGV" 
SCRIPTS="$ROOT/scripts"

python "$SCRIPTS/generate_igv_sessions.py" \
    --peaks "$ROOT/figures/fig1/inputs/123_peaks_annotated-bobsent-20260825.xlsx" \
    --out-dir "$IGVDIR" \
    --top-n 20 \
    --min-reps 2
```

output: 

```
Top 20 peaks:
Rank    adj_sf  Family                     Accession
   1    0.9858  Tombusviridae              NC_055564.1|2620-2849_intergenic
   2    0.9242  Flaviviridae               NC_077001.1|12448-13220_inferred_3UTR
   3    0.9081  Flaviviridae               NC_015843.2|10322-10990_inferred_3UTR
   4    0.9078  Flaviviridae               NC_034017.1|10335-10884_3UTR
   5    0.9077  Solemoviridae              NC_021564.1|4129-4350_intergenic
   6    0.8981  Flaviviridae               NC_034151.1|10381-10937_3UTR
   7    0.8851  Botourmiaviridae           NC_076379.1|2658-2844_inferred_3UTR
   8    0.8831  Botourmiaviridae           NC_076756.1|1992-2420_inferred_3UTR
   9    0.8764  Hypoviridae                NC_076662.1|13267-14211_inferred_3UTR
  10    0.8750  Tornidovirineae            NC_046963.1|27979-28189_intergenic
  11    0.8744  Flaviviridae               NC_008719.1|10284-10793_inferred_3UTR
  12    0.8735  Flaviviridae               NC_012534.1|10325-10941_inferred_3UTR
  13    0.8725  Flaviviridae               NC_038430.1|9355-9609_3UTR
  14    0.8697  Secoviridae                NC_078582.1|5640-7083_3UTR
  15    0.8688  Flaviviridae               NC_025677.1|12480-12983_3UTR
  16    0.8683  Flaviviridae               NC_012932.1|10072-11064_inferred_3UTR
  17    0.8681  Flaviviridae               NC_012533.1|10283-10723_inferred_3UTR
  18    0.8650  Solemoviridae              NC_075424.1|5865-6227_inferred_3UTR
  19    0.8636  Flaviviridae               NC_005039.1|10378-10857_inferred_3UTR
  20    0.8634  Solemoviridae              NC_034836.1|5848-6176_inferred_3UTR

Run from project root:
  bash /Users/jeaninegezelle/Desktop/viromelib/figures/IGV/run_igv_sessions.sh
```