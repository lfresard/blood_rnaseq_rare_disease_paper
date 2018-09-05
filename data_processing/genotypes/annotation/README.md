# Annotate and format vcfs files

## Getting Started

These instructions will allow to annotate vcf files for gnomAD AF and CADD score and to process them for outlier analysis

Files you need to get and transform

## 1- create a vcf file containg both SNPs and Indels CADD annotation

```
zcat cadd_v1.3.vcf.gz InDels.v1.2.noheader.vcf.gz| bgzip >CADD_SNP.INDELS.vcf.gz
vcf-sort -c CADD_SNP.INDELS.vcf.gz | bgzip > CADD_SNP.INDELS.sorted.vcf.gz
tabix CADD_SNP.INDELS.sorted.vcf.gz
```

## 2- Transform vcf files in homogenized format for further analysis
```
find /users/xli6/projects/udn/variant_filter/vcf_processed/filter_byind/*.filtered.vcf.gz  | parallel  "basename {} .filtered.vcf.gz" | parallel  -j 8 "bash /users/lfresard/repos/rare_disease/scripts/vcf_handling/homogenize_vcf_single_sample.sh /users/xli6/projects/udn/variant_filter/vcf_processed/filter_byind/{}.filtered.vcf.gz {} /srv/scratch/restricted/rare_diseases/data/vcfs/for_freeze_exac_filt" > log_normalize_freeze_exac_pipeline_2018_06_13.txt 2>&1 &
```

## 3- Annotate for CADD score and gnomAD AF
```
bash /users/lfresard/repos/rare_disease/scripts/vcf_handling/launch_vcf_anno_chr.sh > log_vcf_gnomad_anno_2018_06_18.txt 2>&1 &
```

## 4- Filter vcf for Rare Variants and output in bed file
* Filter for variants with Allele frequency <= 0.01 (keeps singletons)
* Bedtools intersect with gene bed file to get gene name associated with rare variant

```
bash /users/lfresard/repos/rare_disease/scripts/vcf_handling/launch_RD_vcf_to_RVbed.sh /srv/scratch/restricted/rare_diseases/data/vcfs/for_freeze >log_RD_vcf_to_RVbed_nf_2018_02_16.txt 2>&1 &
```
## 5- Filter for rare variants near annotated splicing junctions
This data is used in figure 3c to estimate the number of rare variants that are within 20bp of splicing junctions

* bed_files_dir=/srv/scratch/restricted/rare_diseases/data/bed/for_freeze/
* RV_directory=/srv/scratch/restricted/rare_diseases/analysis/outlier_analysis/splicing/for_freeze
* sample_script=/users/lfresard/repos/rare_disease/scripts/splicing_analysis/outlier/process_sample_RV_annotated_junctions_window.sh

```
bash /users/lfresard/repos/rare_disease/scripts/splicing_analysis/outlier/RV_near_junction_window.sh > log_RV_nearjunction_2018_04_25.txt 2>&1
```
