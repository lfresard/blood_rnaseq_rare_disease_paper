# Process new batch of RNA-seq data for rare disease project

From fastq files to gene expression data

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Steps
* Adaptor trimming (cutadapt)
* Mapping (STAR)
* Filter for ampping quality >=30 and remove PCR duplicates (Picard)
* Quantify

Take batch number as variable.
Need to create output directories before running script
FASTQ_DIR=/srv/scratch/restricted/rare_diseases/data/fastq/batch${batch_number}
BAM_DIR=/srv/scratch/restricted/rare_diseases/data/mapping/batch${batch_number}
EXP_DIR=/srv/scratch/restricted/rare_diseases/data/quantification/rsem

```
bash /users/lfresard/repos/rare_disease/scripts/fastq_handling/RD_analysis.sh 9 > log_RDbatch9_analysis_2018_03_22.txt 2>&1 &
```
