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
find $filtered_variant_dir/*.filtered.vcf.gz  | parallel  "basename {} .filtered.vcf.gz" | parallel  -j 8 "bash homogenize_vcf_single_sample.sh $filtered_variant_dir/{}.filtered.vcf.gz {} $homogenized_vcf_dir" > log_file.txt 2>&1 &
```

## 3- Annotate for CADD score and gnomAD AF
```
bash launch_vcf_anno_chr.sh > log_file_name.txt 2>&1 &
```

## 4- Filter vcf for Rare Variants and output in bed file
* Filter for variants with Allele frequency <= 0.01 (keeps singletons)
* Bedtools intersect with gene bed file to get gene name associated with rare variant

```
bash launch_RD_vcf_to_RVbed.sh $annotated_vcf_dir >log_file_name.txt 2>&1 &
```
## 5- Filter for rare variants near annotated splicing junctions
This data is used in figure 3c to estimate the number of rare variants that are within 20bp of splicing junctions

```
bash RV_near_junction_window.sh > log_file_name.txt 2>&1
```
