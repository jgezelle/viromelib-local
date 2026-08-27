```bash
conda activate viromelib

ROOT=~/Desktop/viromelib
cd $ROOT
FA=$ROOT/data/pilot_library_20260319.uniqueIDs.fa
CMDIR=$ROOT/work/06_infernal/cm
TBLOUT=$ROOT/work/06_infernal/tblout
mkdir -p $TBLOUT

CMS=(
    "CL_Class_1a_xrRNA_final.cm"
    "CL_Class_1b_xrRNA_final.cm"
    "CL_Class_2_xrRNA_final.cm"
)

for CM in "${CMS[@]}"; do
    echo "Searching: $CM"
    cmsearch --tblout "$TBLOUT/${CM}.tbl" \
             --cpu 8 \
             "$CMDIR/$CM" "$FA"
    echo "Done: $CM"
done
```

Try the same, but for flaviLib: 

```bash
conda activate viromelib

ROOT=~/Desktop/viromelib
cd $ROOT
FA=$ROOT/data/flavi_library_20260228.uniqueIDs.fa
CMDIR=$ROOT/work/06_infernal/cm
TBLOUT=$ROOT/work/06_infernal/flaviLib/tblout
mkdir -p $TBLOUT

CMS=(
    "CL_Class_1a_xrRNA_final.cm"
    "CL_Class_1b_xrRNA_final.cm"
    "CL_Class_2_xrRNA_final.cm"
)

for CM in "${CMS[@]}"; do
    echo "Searching: $CM"
    cmsearch --tblout "$TBLOUT/${CM}.tbl" \
             --cpu 8 \
             "$CMDIR/$CM" "$FA"
    echo "Done: $CM"
done
```

take into `viromelib/work/cm_distributions_viromeLib_flaviviridae.Rmd`
for distribution analysis with proper file paths 
