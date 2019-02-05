# Process new batch of RNA-seq data for rare disease project

From fastq files to gene expression data

### Steps
* Adaptor trimming (cutadapt)
* Mapping (STAR)
* Filter for mapping quality >=30 and remove PCR duplicates (Picard)
* Quantify

Take batch number as variable.
Need to create output directories before running script
FASTQ_DIR=$<fastq_dir>/batch${batch_number} .
BAM_DIR=<$bam_dir>/batch${batch_number} . 
EXP_DIR=<$expression_dir>/rsem . 

```
bash RD_analysis.sh <batch_number> > log_file_name.txt 2>&1 &
```
