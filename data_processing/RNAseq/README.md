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
FASTQ_DIR=$<fastq_dir>/batch${batch_number}\\
BAM_DIR=<$bam_dir>/batch${batch_number}\\
EXP_DIR=<$expression_dir>/rsem\\

```
bash RD_analysis.sh <batch_number> > log_file_name.txt 2>&1 &
```
