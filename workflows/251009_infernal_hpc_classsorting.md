# infernal searches 
## against the Xrn1-stop peak windows .fa, using established xrRNA CMs 
### **organize workspace:**
```bash
ROOT=/groups/as6282_gp/scratch_bkup/jgg2144/xrrna_seq/viromelib
# original location of the peak windows .fa output 
FA_DIR="$ROOT/work/04_align/fullacc" && cd $FA_DIR
WORK="$ROOT/work/06_infernal"
mkdir -p "$WORK/tblout" "$WORK/cm" "$WORK/fa" "$WORK/groups_out" # where outputs and inputs will go

# either copy peak windows output to infernal dir or give path
cp 250904-r0368-X-01-A3405_norm_PASS_peak_windows.fa ../../06_infernal/fa/ && cd $WORK
```

```bash
# if not in conda env 
conda activate viromelib 
```

### **Run infernal search against hits for all known cms:**
```bash
# within 06_infernal:
cmsearch --tblout tblout/cl4_r368_norm.tbl cm/cl4.cm fa/250904-r0368-X-01-A3405_norm_PASS_peak_windows.fa

cmsearch --tblout tblout/cl3_r368_norm.tbl cm/cl3.cm fa/250904-r0368-X-01-A3405_norm_PASS_peak_windows.fa

cmsearch --tblout tblout/cl2_r368_norm.tbl cm/cl2.cm fa/250904-r0368-X-01-A3405_norm_PASS_peak_windows.fa

cmsearch --tblout tblout/cl1b_r368_norm.tbl cm/cl1b.cm fa/250904-r0368-X-01-A3405_norm_PASS_peak_windows.fa

cmsearch --tblout tblout/cl1a_r368_norm.tbl cm/cl1a.cm fa/250904-r0368-X-01-A3405_norm_PASS_peak_windows.fa
```
## Build run map summary csv 
Python script: <br>
**buildInfernalRunMap.py** <br>
script to create a single mapping CSV where each row corresponds to one cmsearch --tblout file. <br>
*infernal_runs.csv* **can then be used for the cm class sorting script sortclasses.py**

Each row contains:
- tbl_path:   path to a tblout file
- dataset:    chosen dataset label (e.g., r0368, r0366_combined)
- fasta_path: path to the FASTA used to define window IDs for this dataset (e.g., peak_windows.fa)
- query_name: the query token found inside the tbl (usually the CM name printed by cmsearch)
- class_label: a canonical label among {class1a, class1b, class2, class3, class4} (auto-filled if possible or left blank)
- map_source   : 'explicit' (from map), 'heuristic', or 'blank'

- For each dataset block given the script writes one CSV row per tbl output from infernal and fills dataset and fasta_path,
- then sets class_label using a simple map provided via --cm-map 

Notes:
- The query is taken from column 3 (index 2) of the first non-comment data line.
- If column 3 is missing or '-', the parser looks ahead in the first 6 tokens for a likely query.
- Mapping order: explicit --cm-map pairs (case-insensitive key match) -> optional CSV map -> optional heuristic.

### Example run
--------
Two datasets, 5 tbls each. <br>
```py
SCRIPT="$ROOT"/scripts/buildInfernalRunMap.py
# make sure upstream files are uniquely named from other runs 
# cm-map where to enter the query column name from .tbl
python "$SCRIPT" \
  --add-dataset r0368 fa/250904-r0368-X-01-A3405_norm_PASS_peak_windows.fa "tblout/*_r368_norm.tbl" \
  --cm-map cl4=class4 \
  --cm-map class3=class3 \
  --cm-map cl2=class2 \
  --cm-map class1bresult=class1b \
  --cm-map class1a_result=class1a \
  --out csv/251009_infernal_runs.csv # add date YYMMDD
  ```

### Example output 
`[OK] wrote /groups/as6282_gp/scratch_bkup/jgg2144/xrrna_seq/viromelib/work/06_infernal/251009_infernal_runs.csv with 5 rows
[INFO] mapping sources — explicit: 5, heuristic: 0`

## Determine class (cm) distribution among hits: 
Python script: <br>
**sortclasses.py** <br>
Window-level summarizer for Infernal cmsearch --tblout results, using a single mapping CSV.

**Input** CSV (one row per tbl): <br>
  `tbl_path,dataset,fasta_path,query_name,class_label` <br>
  - tbl_path    : path to the .tbl file (should be dir `tblout`)
  - dataset     : dataset label (e.g., `r0368`, `r0366_combined`)
  - fasta_path  : path to FASTA listing the windows for that dataset (should be dir `fa`)
  - query_name  : query token found in the tbl (for reference)
  - class_label : one of `{class1a,class1b,class2,class3,class4}`

**Outputs** (by dataset): out dir should be `groups_out` <br>
  `<out>/<dataset>/window_level/class_<class>.txt`; can change name to just `<class>.txt` too but have to change in next script as well <br>
  `<out>/<dataset>/window_level/multiclass.csv` <br>
  `<out>/<dataset>/window_level/no_group.txt` <br>

**Outputs** (across datasets): out dir should be `groups_out` <br>
  `<out>/window_hits_long.csv`
  `<out>/window_counts_by_class.csv`
  `<out>/summary.csv`

Notes
-----
- peak-window level only. Keeps multiple peak windows per accession, based on the passing peaks.
    - **filtering is an upstream step**
    - **runs .csv is an upstream step**
- "Unique" means removing exact duplicate rows so a window is counted once per class per dataset.
- Optional i-Evalue filter (best-effort from numeric tokens near the end of each line); acts only if a value is found.

### Example Run:
```bash
  python "$ROOT"/scripts/sortclasses.py \
    --runs-csv 251009_infernal_runs.csv \
    --out-dir groups_out \
    --i-evalue 1e-2
```
### Example output: 
`[OK] wrote outputs under: /groups/as6282_gp/scratch_bkup/jgg2144/xrrna_seq/viromelib/work/06_infernal/groups_out` <br>
`Per-dataset outputs are under <out>/<dataset>/window_level/`

## Merge cm class info with the PASS.csv following filtering, to consolidate all information into one csv
Python script: <br>
**joinPASScsv_Infernal.py** <br>

Make one sheet that merges:
  - class membership (from class text files; this follows the infernal output and **sortclasses.py workflow**)
  - window sequence (from the peak windows FASTA)
  - PASS metrics from PASS CSV (post-filter), matched by columns: accession, peak_position

**Inputs:** <br>
`--pass-csv`         : PASS table for this dataset (must contain 'accession' and either 'peak_position'
                     or 'peak_positions' [comma-separated ints]) <br>
`--fasta`            : peak windows FASTA for this dataset (headers are the window IDs) <br>
`--class-file`       : Repeatable. Format `CLASS=PATH` (e.g., class1a=class_class1a.txt); Use 'no_group' as a class name to include no_group.txt.<br>
`--dataset`          : dataset label to record in the output (e.g., r0368) <br>
`--out-csv`          : output CSV path <br>

Output columns (in this order):
  `dataset, class_label, window_id, accession, peak, pos, len, sequence,`  [then `PASS` columns if present]


**Example run:** <br>
```bash
cd groups_out

python "$ROOT"/scripts/joinPASScsv_Infernal.py \
  --pass-csv "$ROOT/work/04_align/fullacc/250904-r0368-X-01-A3405_norm_filtered.PASS.csv" \
  --fasta    ../fa/250904-r0368-X-01-A3405_norm_PASS_peak_windows.fa \
  --class-file class1a=r0368/window_level/class_class1a.txt \
  --class-file class1b=r0368/window_level/class_class1b.txt \
  --class-file class2=r0368/window_level/class_class2.txt \
  --class-file class3=r0368/window_level/class_class3.txt \
  --class-file class4=r0368/window_level/class_class4.txt \
  --class-file no_group=r0368/window_level/no_group.txt \
  --dataset  r0368 \
  --out-csv  251009_r0368_infernal_PASS_join.csv # add date YYMMDD to name 
  ```
## scp to locally analyze cm class distribution
```bash
# on local machine in viromelib dir:
scp -r jgg2144@hpc.c2b2.columbia.edu:/groups/as6282_gp/scratch_bkup/jgg2144/xrrna_seq/viromelib/work/06_infernal/groups_out .
```