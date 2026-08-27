fetch info from NCBI 
login with ORCID , jeaninegezelle@gmail.com, UN jgezelle, go to account settings, scroll to bottom, get API key: `938f48a619230ed0fca54316f3ac275e8508`

```bash 
# within conda env viromelib, local 
python fetch_ncbi_family_counts.py virus_stats.csv ncbi_family_counts.csv --api-key 938f48a619230ed0fca54316f3ac275e8508

```
Now set NCBI_COUNTS_FILE (in the R script) and re-run:
```bash
NCBI_COUNTS_FILE <- "ncbi_family_counts.csv"

in R notebook:
```
source("virome_fig1a_pipeline.R")
```