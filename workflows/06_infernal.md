## Jul 20, 2026
build class 4 cm 
reference - Dilweg et al. 2023 Sci Rep "Xrn1-resistant RNA motifs are disseminated throughout the RNA virome and are able to block scanning ribosomes"

saved sequences as fasta under `fa` folder. 

```bash
conda activate viromelib
ROOT=~/Desktop/viromelib/work/06_infernal
cd $ROOT 
mlocarna --stockholm --consensus-structure=alifold fa/20260720_cl4_input_mlocarna.fa
mv fa/20260720_cl4_input_mlocarna.out stockholm/20260720_cl4_input_mlocarna.out
# rename output to contain cl4 in name 
```

infernal 
```bash
cmbuild cm/cl4_20260720-2.cm stockholm/cl4_20260720_result.stk # manually added base pair in second stem. 
cmcalibrate cm/cl4_20260720-2.cm
cmsearch --tblout tblout/cl4_20260720-2.tbl cm/cl4_20260720-2.cm fa/RNAvirus_CompleteAssem_No_CoV-2.fa
```
### Beny seqs have an e-value greater than 1... need to improve CM 

## `cmalign` 
in dir `06_infernal`: 
```bash
mkdir -p cmalign
grep -v "^#" tblout/cl4_20260720-2.tbl | grep -v "^$" | \
  awk '{print $1"/"$8"-"$9, $8, $9, $1}' > cmalign/hits.fetchlist

# actually, do awk only for positive strand (col 10) and greater than 2 e-value (col 17)
grep -v "^#" tblout/cl4_20260720-2.tbl | grep -v "^$" | \
  awk '$10=="+" && $16<=2 {print $1"/"$8"-"$9, $8, $9, $1}' > cmalign/hits.fetchlist
```
need to install esl-fetch.. for some reason not included in my infernal install? use hmmer
```bash
conda install -n viromelib -c bioconda hmmer # done once 
```

now, run `esl-sfetch` which will generate the new .fa 
```bash
esl-sfetch --index fa/RNAvirus_CompleteAssem_No_CoV-2.fa # index, do once for search fa 

esl-sfetch -Cf fa/RNAvirus_CompleteAssem_No_CoV-2.fa cmalign/hits.fetchlist > cmalign/hits.fa
# get counts
grep -c "^>" cmalign/hits.fa
```
17 hits 

### run `cmalign`
```
cmalign --cpu 4 -o cmalign/hits_aligned.stk cm/cl4_20260720-2.cm cmalign/hits.fa
cmbuild -F -n cl4_v3 cm/cl4_v3.cm cmalign/hits_aligned.stk
cmcalibrate --cpu 4 cm/cl4_v3.cm
cmsearch --tblout tblout/cl4_v3.tbl --cpu 4 cm/cl4_v3.cm fa/RNAvirus_CompleteAssem_No_CoV-2.fa
```

42 seqs, got slightly better. 
add some diversity with the 'xrRNAs everywhere' 
```bash
cat > cmalign/new_accessions.txt << 'EOF'
NC_004750
NC_003516
KX665538
EU480547
KP316705
KP316716
KP316688
EU480523
DQ462122
KP316974
NC_003513
KP316967
KP316935
KP316952
KP006484
KP316958
KP316957
KP316937
AY700056
NC_039224.1
NC_039227.1
NC_040760.1
NC_003837.1
NC_003836.1
NC_003838.1
NC_003811.1
KM434206
KX013371.1
KX013372.1
MF170158.1
AY775057.1
MF170160.1
NC_002040.1
LC380680.1
NC_012136.1
KU955493
MF061245.1
JN127346
NC_018448.1
MG995734.1
MH036372.1
NC_023016.1
NC_023017.1
NC_039002.1
NC_024808
MG646923.1
NC_003603.1
NC_029783.1
EU099844
NC_003510
EOF

efetch -db nuccore -input cmalign/new_accessions.txt -format fasta > cmalign/new_seqs.fa
grep -c "^>" cmalign/new_seqs.fa
```

```bash
# new cmsearch with those accessions
cmsearch --tblout cmalign/new_seqs.tbl --cpu 4 cm/cl4_v3.cm cmalign/new_seqs.fa

# sort by score .. 
grep -v "^#" cmalign/new_seqs.tbl | grep -v "^$" | \
  awk '$10=="+"' | \
  sort -k16 -rn | \
  awk '{printf "%-20s %6s %6s  score:%-7s  E:%-10s  %s\n", $1,$8,$9,$16,$17,$18}' | head -40
```

be conservative, get ! (sig) hits only; 
```bash
grep -v "^#" cmalign/new_seqs.tbl | grep -v "^$" | \
  awk '$10=="+" && $17=="!" {print $1, $8, $9, $15, $16}' | \
  sort -k1,1 -k5,5n | \
  awk '!seen[$1]++ {print $1"/"$2"-"$3, $2, $3, $1}' \
  > cmalign/new_hits.fetchlist

wc -l cmalign/new_hits.fetchlist
```
next search: 
```bash
esl-sfetch --index cmalign/new_seqs.fa
esl-sfetch -Cf cmalign/new_seqs.fa cmalign/new_hits.fetchlist > cmalign/new_hits.fa
grep -c "^>" cmalign/new_hits.fa

cat cmalign/hits.fa cmalign/new_hits.fa > cmalign/combined_hits.fa

cmalign --cpu 4 -o cmalign/combined_aligned.stk cm/cl4_v3.cm cmalign/combined_hits.fa

cmbuild -F -n cl4_v4 cm/cl4_v4.cm cmalign/combined_aligned.stk
cmcalibrate --cpu 4 cm/cl4_v4.cm
cmsearch --tblout tblout/cl4_v4.tbl --cpu 4 cm/cl4_v4.cm fa/RNAvirus_CompleteAssem_No_CoV-2.fa
```
move forward with the ! hits 
```bash
grep -v "^#" tblout/cl4_v4.tbl | grep -v "^$" | \
  awk '$10=="+" && $17=="!" {print $1"/"$8"-"$9, $8, $9, $1}' > cmalign/v4_sig_hits.fetchlist

wc -l cmalign/v4_sig_hits.fetchlist
```
make 5th version of cm-- 
```bash
esl-sfetch -Cf fa/RNAvirus_CompleteAssem_No_CoV-2.fa cmalign/v4_sig_hits.fetchlist > cmalign/v4_sig_hits.fa
grep -c "^>" cmalign/v4_sig_hits.fa
#28 
# new hits in fa to be incorporated by cmalign
cmalign --cpu 4 -o cmalign/v4_sig_aligned.stk cm/cl4_v4.cm cmalign/v4_sig_hits.fa

cmbuild -F -n cl4_v5 cm/cl4_v5.cm cmalign/v4_sig_aligned.stk

cmcalibrate --cpu 4 cm/cl4_v5.cm

cmsearch --tblout tblout/cl4_v5.tbl --cpu 4 cm/cl4_v5.cm fa/RNAvirus_CompleteAssem_No_CoV-2.fa
```
do one more iteration with just the 18 `!` hits 
```bash
grep -v "^#" tblout/cl4_v5.tbl | grep -v "^$" | \
  awk '$10=="+" && $17=="!" {print $1"/"$8"-"$9, $8, $9, $1}' > cmalign/v5_sig_hits.fetchlist

wc -l cmalign/v5_sig_hits.fetchlist # expect 18 with the results 
```
```bash 
esl-sfetch -Cf fa/RNAvirus_CompleteAssem_No_CoV-2.fa cmalign/v5_sig_hits.fetchlist > cmalign/v5_sig_hits.fa
grep -c "^>" cmalign/v5_sig_hits.fa
```
```bash
cmalign --cpu 4 -o cmalign/v5_sig_aligned.stk cm/cl4_v5.cm cmalign/v5_sig_hits.fa

cmbuild -F -n cl4_v6 cm/cl4_v6.cm cmalign/v5_sig_aligned.stk

cmcalibrate --cpu 4 cm/cl4_v6.cm

cmsearch --tblout tblout/cl4_v6.tbl --cpu 4 cm/cl4_v6.cm fa/RNAvirus_CompleteAssem_No_CoV-2.fa
```
one more iteration, even though model seems to have converged on same hits. 

```bash
grep -v "^#" tblout/cl4_v6.tbl | grep -v "^$" | \
  awk '$10=="+" && $17=="!" {print $1"/"$8"-"$9, $8, $9, $1}' > cmalign/v6_sig_hits.fetchlist

wc -l cmalign/v6_sig_hits.fetchlist

esl-sfetch -Cf fa/RNAvirus_CompleteAssem_No_CoV-2.fa cmalign/v6_sig_hits.fetchlist > cmalign/v6_sig_hits.fa

cmalign --cpu 4 -o cmalign/v6_sig_aligned.stk cm/cl4_v6.cm cmalign/v6_sig_hits.fa

cmbuild -F -n cl4_v7 cm/cl4_v7.cm cmalign/v6_sig_aligned.stk

cmcalibrate --cpu 4 cm/cl4_v7.cm

cmsearch --tblout tblout/cl4_v7.tbl --cpu 4 cm/cl4_v7.cm fa/RNAvirus_CompleteAssem_No_CoV-2.fa
```
----------
note, currently, `UUAAACUAC` motif seems to dominate the search. 
could use it as a pre-filter in the RNA virus database, then build the CM only on the flanking structure.. 
search for the DNA version

----------

### align the last iteration 

```bash 
grep -v "^#" tblout/cl4_v7.tbl | grep -v "^$" | \
  awk '$10=="+" && $17=="!" {print $1"/"$8"-"$9, $8, $9, $1}' > cmalign/v7_sig_hits.fetchlist

esl-sfetch -Cf fa/RNAvirus_CompleteAssem_No_CoV-2.fa cmalign/v7_sig_hits.fetchlist > cmalign/v7_sig_hits.fa

cmalign --cpu 4 -o cmalign/v7_sig_aligned.stk cm/cl4_v7.cm cmalign/v7_sig_hits.fa
```

### tried the cm search against peak windows from virome lib hits secoviridae family only  - no hits 
```bash
 cmsearch --tblout tblout/cl4_v7_virome-seco-peaks.tbl --cpu 4 \
  cm/cl4_v7.cm fa/Secoviridae.peak_windows.merge.fa

```
### webserver R-scape results 
input `v7_sig_aligned.stk`

try two runs: 1) Evaluate region for conserved structure, 2) Improve given structure with CaCoFold 

CaCoFold made stem two four bases longer. 

5 sig. covarying bp 

### install R-scape
#### manually click download at `https://github.com/EddyRivasLab/R-scape/releases/tag/v2.6.11`

```bash
cd ~/Downloads
tar xzf rscape_v2.6.11.tar.gz
cd rscape_v2.6.11

./configure --prefix=$HOME/local
make -j 4
make install

echo 'export PATH=$HOME/local/bin:$PATH' >> ~/.zshrc
source ~/.zshrc

which R-scape
```

### use R-scape locally:
```bash
cd ~/Desktop/viromelib/work/06_infernal
conda activate viromelib

mkdir -p rscape/cl4_v7
# mode: 
# From the manual, "Improve given structure with CaCoFold" corresponds to the two-set test + cacofold combination:
R-scape -s --cacofold \
  --outdir rscape/cl4_v7 \
  --outname cl4_v7_improve \
  cmalign/v7_sig_aligned.stk

```

From the manual:

`-s` alone = "Evaluate given structure" (two-set test against existing SS_cons)

`--cacofold` alone = "Predict new structure" (ignores SS_cons)

`-s --cacofold` together = "Improve given structure with CaCoFold" 

(two-set test + uses SS_cons as starting point + CaCoFold improves it)

The web server "Improve" mode runs both statistical tests (one on the annotated pairs in your SS_cons, one on all other possible pairs), then CaCoFold builds the best structure incorporating all significant pairs from both tests. 

### improve cm with the cacofold/R-scape output:
```bash
# check the new SS_cons
grep "SS_cons" rscape/cl4_v7/cl4_v7_improve.cacofold.sto

# build new CM from it using --hand to respect the annotated structure
cmbuild --hand -F cm/cl4_v8-rscape.cm rscape/cl4_v7/cl4_v7_improve.cacofold.sto

cmcalibrate --cpu 4 cm/cl4_v8-rscape.cm

cmsearch --tblout tblout/cl4_v8-rscape.tbl --cpu 4 cm/cl4_v8-rscape.cm fa/RNAvirus_CompleteAssem_No_CoV-2.fa
```

#### get fetchlist with new `!` hits: 
```bash
grep -v "^#" tblout/cl4_v8-rscape.tbl | grep -v "^$" | \
  awk '$10=="+" && $17=="!" {print $1"/"$8"-"$9, $8, $9, $1}' > cmalign/v8_sig_hits.fetchlist

wc -l cmalign/v8_sig_hits.fetchlist # 29 sig hits 
```
#### run `cmalign`: 
```bash
esl-sfetch -Cf fa/RNAvirus_CompleteAssem_No_CoV-2.fa cmalign/v8_sig_hits.fetchlist > cmalign/v8_sig_hits.fa

grep -c "^>" cmalign/v8_sig_hits.fa

cmalign --cpu 4 -o cmalign/v8_sig_aligned.stk cm/cl4_v8-rscape.cm cmalign/v8_sig_hits.fa

cmbuild --hand -F -n cl4_v9 cm/cl4_v9.cm cmalign/v8_sig_aligned.stk

cmcalibrate --cpu 4 cm/cl4_v9.cm

cmsearch --tblout tblout/cl4_v9.tbl --cpu 4 cm/cl4_v9.cm fa/RNAvirus_CompleteAssem_No_CoV-2.fa
```
### the model has converged to 30 `!` hits , stop iterating 
### Run `R-scape` again
```bash
mkdir -p rscape/cl4_v9
R-scape -s --cacofold \
  --outdir rscape/cl4_v9 \
  --outname cl4_v9_improve \
  cmalign/v8_sig_aligned.stk
```


save all to github
```bash
cd ~/Desktop/viromelib
git init
```
need to update `.gitignore`: 

also, `flaviLib` dir `$ROOT/flaviLib` has its own git repo already, so I will point to it as a submodule. 

```bash
cd ~/Desktop/viromelib/flaviLib
```
```bash
cd ~/Desktop/viromelib
git rm --cached -f flaviLib

git submodule add git@github.com:jgezelle/flaviLib.git flaviLib

cat .gitmodules  # should show the submodule entry
```
```
[submodule "flaviLib"]
        path = flaviLib
        url = git@github.com:jgezelle/flaviLib.git
```

reorganize directory then git add: 

```bash
# remove 
git rm -f -r --cached .

# add specific files i want 
git add .gitignore .gitmodules flaviLib
git add work/06_infernal/cm/
git add work/06_infernal/06_infernal.md
git add work/*.md
git add work/*.Rmd
git add work/*.rmd
git add htpmut/

# check 
git status --short | grep "^A"
```

### finally ready to commit: 
```bash
git commit -m "initial commit with infernal pipeline, CM files, analysis notebooks"

# create the repo on github.com first if not done yet, then:
git remote add origin git@github.com:jgezelle/viromelib-local.git

git branch -M main
git push -u origin main

```

**move this .md into the `workflows` subdir of `viromelib`**