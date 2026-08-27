Initialize
```bash
replicates=("1" "2" "3")
peak_script="scripts/filterPeaks_normtoCTRL_v7.3.py"
calibrate_script="scripts/calibratePeakNull.py"
outdir="work/06_filter_v7"
fragment_dir="work/05_fragments"
tiles_file="work/05_fragments/pilot_library_20260319.uniqueIDs.bed"
```

Create fragment files

```bash
for f in /groups/as6282_gp/scratch_bkup/jgg2144/xrrna_seq/viromelib/work/04_align/fs4/*.bam; do
    echo $f
    python3 scripts/make_fragment_file.py $f > work/05_fragments/$(basename $f .bam)
done
```

Calibrate

```bash
for r in "${replicates[@]}"; do 
    python3 $calibrate_script $fragment_dir/*-X${r}*k1.mq0.csv $tiles_file $outdir/${r}_calibrate '|' tee $outdir/${r}_calibrate.out.txt
done
```

From this we get the three alpha parameter point estimates for each of the replicates.

```bash
alphas=("2.0885" "1.7196" "2.1398")


# for r in $replicates; do 
for i in ${!replicates[@]}; do
    r=${replicates[$i]}
    alpha=${alphas[$i]}
    echo "Processing replicate: $r"
    python $peak_script \
        $fragment_dir/*X${r}-*k1.mq0.csv $fragment_dir/*R${r}-*k1.mq0.csv $tiles_file $outdir/${r}_peaks.csv  \
        --out-accession-csv $outdir/${r}_accessions.csv \
        --out-gtf $outdir/${r}_peaks.gtf \
        --null negbinom \
        --nb-alpha $alpha
done
```

Perform CM hit analysis

```bash
cms=(
    "CL_Class_1a_xrRNA_final.cm"
    "CL_Class_1b_xrRNA_final.cm"
    "CL_Class_2_xrRNA_final.cm"
    "cl3_NAR_2025.align.mpkadded.edited.confident.dedup.IDcurated.SS_cons.cm"
    "cl4.cm"
)
cmdir="data/jeanine_cm"
ref=data/pilot_library_20260319.uniqueIDs.fa
tblout=work/08_infernal
for cm in "${cms[@]}"; do
    echo "Processing CM: $cm"
    cmsearch --tblout $tblout/$cm.tbl $cmdir/$cm $ref
done
```

Combine replicate results

```bash
replicate_csvs=$(printf -- "$outdir/%s_peaks.csv " "${replicates[@]}")
python3 scripts/combine_replicates.py $replicate_csvs --labels R1 R2 R3 --out $outdir/123_peaks_combined.csv
```

Annotate peak files

```bash
annout=work/10_postproc
mkdir -p $annout
infernal_args=$(printf -- "--infernal-tbl $tblout/%s.tbl " "${cms[@]}")
# maps to --infernal-tbl data/jeanine_cm/CL_Class_1a_xrRNA_final.cm --infernal-tbl data/jea...virus
for r in ${replicates[@]}; do 
    echo "Annotating peaks for replicate: $r"
    python scripts/post_process_peaks.py $outdir/${r}_peaks.csv \
        --out-raw $annout/${r}_peaks_annotated.csv \
        --out-collapsed $annout/${r}_peaks_per_accession.csv \
        --min-treated-sf 0.15 \
        --min-reads 20 \
        --family-csv data/virus_stats.csv \
        $infernal_args
done

python scripts/post_process_peaks.py $outdir/123_peaks_combined.csv \
        --out-raw $annout/123_peaks_annotated.csv \
        --out-collapsed $annout/123_peaks_per_accession.csv \
        --family-csv data/virus_stats.csv \
        $infernal_args

python3 scripts/extract_downstream_seqs.py $outdir/123_peaks_combined.csv --fasta  data/pilot_library_20260319.uniqueIDs.fa  --out $outdir/123_peaks_combined_downstream_seqs_unique.fa
```


Make bedpe files

```bash
tracks_dir=work/11_tracks
mkdir -p $tracks_dir
for r in ${replicates[@]}; do
    python scripts/fragments_to_tracks.py $fragment_dir/*X$r-*k1.mq0.csv \
        --tile-bed $tiles_file \
        --out-bedgraph $tracks_dir/treated_R$r.5p.bedgraph \
        --out-tile-bed $tracks_dir/treated_R$r.tiles.bed \
        --track-name "treated_R$r"
    python scripts/fragments_to_tracks.py $fragment_dir/*R$r-*k1.mq0.csv \
        --tile-bed $tiles_file \
        --out-bedgraph $tracks_dir/control_R$r.5p.bedgraph \
        --out-tile-bed $tracks_dir/control_R$r.tiles.bed \
        --track-name "control_R$r"
done
```

Make an IGV session:

```bash
python scripts/make_igv_session.py \
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
    --gtf work/06_filter_v7/{1,2,3}_peaks.gtf \
    --bed work/05_fragments/pilot_library_20260319.uniqueIDs.bed \
    --locus "NC_024114.1|2493-2824_inferred_3UTR:0-196" \
    --out 3replicates_generated.xml


```

Run locarna
```bash
mkdir $outdir/seq_split
cd $outdir/seq_split
# proprietary util from me, can be replaced with another tool that splits multi-FASTA into single-FASTA files named by the sequence ID.  Output files are named like NC_038434.1|10194-10465_3UTR|191.fa
fasta_split ../123_peaks_combined_downstream_seqs.fa 
mkdir locarna_out
for x in *.fa; do
    for y in *.fa; do
        if [[ "$x" < "$y" ]] || [[ "$x" == "$y" ]]; then
            echo "locarna --stockholm \"locarna_out/${x}--${y}.stk\" --score-components --local-output --sequ-local=1 --struct-local=1 \"$x\" \"$y\" > \"locarna_out/${x}--${y}.out\""
        fi
    done
done > locarna_commands.sh
echo
```

Perform clustering in  scripts/clustering.R

Run mlocarna: 