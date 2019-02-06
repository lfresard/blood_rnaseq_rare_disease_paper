# Rare disease Blood RNA-seq pipeline

This repository contains code to run all components of the Rare disease Blood RNA-seq pipeline in order to highlight candidate genes. 

The following software and data are necessary:

* cutadapt 1.11 https://pypi.python.org/pypi/cutadapt  
* STAR 2.4.0j https://github.com/alexdobin/STAR  
* Python 2.7. https://www.python.org/  
* Picard tools https://broadinstitute.github.io/picard/ and https://github.com/broadinstitute/picard  
* Samtools http://www.htslib.org/   
* RSEM v1.2.21 https://github.com/deweylab/RSEM/   
* [EXAC](http://exac.broadinstitute.org/downloads) data 

Many scripts rely on the existence of a [metadata](./metadata.md) file containing information on samples to analyze.


## Overview  

### 1. Genetic data processing
* 1.1 Variant filtering
* 1.2 VCF Annotation 

### 2. RNA-seq data processing
* 2.1 Adapters trimming 
* 2.2 Alignement
* 2.3 Filters
* 2.4 Quantification
* 2.5 Gene expression data processing 
* 2.6 Splicing data processing 

### 3. Download phenotypic information
* 3.1 Get phenotype to gene link from Human Phenotype Ontology:
* 3.2 Parse HPO database

####  4. Combine information to highlight candidate genes
* 4.1 Expression outlier filtering
* 4.2 Splicing outlier filtering

#### 5. Highlight candidate genes

# 1. Genetic data processing

The variant QC generally follows the ExAC guideline as in [ref1](https://doi.org/10.1016/j.ajhg.2018.05.002) and [ref2](https://doi.org/10.1038/nature19057).

## 1.1 Variant filtering
* keep only "PASS" variants, determined by the GATK variant quality score recalibration (VQSR)

adhoc filter for non-GATK processed vcfs: keep "." or "num prev seen samples > 30"
* as in ["hail"](https://hail.is/docs/stable/index.html) api

```
.filter_genotypes(
         """
         (g.isHomRef && (g.ad[0] / g.dp < 0.8 || g.gq < 20 || g.dp < 20)) ||
         (g.isHomVar && (g.ad[1] / g.dp < 0.8 || g.pl[0] < 20 || g.dp < 20)) ||
         (g.isHet && ( (g.ad[0] + g.ad[1]) / g.dp < 0.8 || g.ad[1] / g.dp < 0.20 || g.pl[0] < 20 || g.dp < 20))
         """, keep = False)
```         


* as in VCF filters:

Filters | FORMAT/GT | FORMAT/DP | FORMAT/GQ | FORMAT/PL | FORMAT/AD
---|---|---|---|---|---
keep|=='0/0'|>20|>20|PL[0]<20|AD[0]/DP>0.8
keep|=='0/1'|>20|>20|PL[1]<20|AD[0]/DP>0.2 & AD[1]/DP>0.2 & (AD[0]+AD[1])/DP>0.2
keep|=='1/1'|>20|>20|PL[2]<20|AD[1]/DP>0.8

note 1: PL = *-10\*log10(likelihood)*, most likely being *PL=0*;  
note 2: FORMAT information is different across sits/samples/institutions, cannot apply uniform filters
* CGS	GT:AD:DP:GQ:PL
* CGS	GT:GQ:PL
* CHEO	GT
* CHEO	GT:AD:DP:GQ:PL
* CHEO	GT:GQ:PL
* CHEO	GT:PL
* CHEO	GT:PL:DP:SP:GQ
* UDN	GT:AD
* UDN	GT:AD:DP
* UDN	GT:AD:DP:GQ:PGT:PID:PL
* UDN	GT:AD:DP:GQ:PL
* UDN	GT:AD:DP:PGT:PID
* UDN	GT:AD:GQ:PGT:PID:PL
* UDN	GT:AD:GQ:PL
* UDN	GT:AD:PGT:PID
* UDN	GT:VR:RR:DP:GQ

Script: [generate_variantfilters.sh](./data_processing/genotypes/variant_filter/generate_variantfilters.sh)

Input: VCF files from each individual

Output: Text file containing list of sites passing filters across individuals

```
bash generate_variantfilters.sh | parallel --jobs 20
```

* Filter vcf files

Input: VCF files and list of sites passing or not filters

Output: Filtered VCF

Script:[filter_allvcf.sh](./data_processing/genotypes/variant_filter/filter_allvcf.sh)

```
bash filter_allvcf.sh | parallel --jobs 20
```

## 1.2 VCF Annotation 
Annotation with [gnomAD] (https://gnomad.broadinstitute.org/downloads) allele frequency  and [CADD score] (https://cadd.gs.washington.edu/download)

### 1.2.1 Transform vcf files in homogenized format for further analysis

Input: Filtered VCF

Output: Filtered VCF in a unified format across samples

Script: [homogenize_vcf_single_sample.sh](./data_processing/genotypes/annotation/homogenize_vcf_single_sample.sh)

```
find $filtered_variant_dir/*.filtered.vcf.gz  | parallel  "basename {} .filtered.vcf.gz" | parallel  -j 8 "bash homogenize_vcf_single_sample.sh $filtered_variant_dir/{}.filtered.vcf.gz {} $homogenized_vcf_dir" > log_file.txt 2>&1 &
```

### 1.2.2 Annotate for CADD score and gnomAD AF
This script will run on 1 chromosome at a time and concatenate results at the end.

Input: Homogenized filtered VCF files from previsous step.

Output: VCF annotated for CADD score and gnomAD allele frequency

Script: [launch_vcf_anno_chr.sh](./data_processing/genotypes/annotation/launch_vcf_anno_chr.sh)
This script uses [vcf_anno_chr.sh](./data_processing/genotypes/annotation/vcf_anno_chr.sh) and [concatenate_gnomadres.sh](./data_processing/genotypes/annotation/concatenate_gnomadres.sh)

```
bash launch_vcf_anno_chr.sh > log_file_name.txt 2>&1 &
```

### 1.2.3 Filter VCF file for Rare Variants and output in bed file
* Filter for variants with Allele frequency <= 0.01 (keeps singletons)
* Bedtools intersect with gene bed file to get gene name associated with rare variant

Script: [launch_RD_vcf_to_RVbed.sh](./data_processing/genotypes/annotation/launch_RD_vcf_to_RVbed.sh)
```
bash launch_RD_vcf_to_RVbed.sh $annotated_vcf_dir >log_file_name.txt 2>&1 &
```
### 1.2.4 Filter for rare variants near annotated splicing junctions
This data is used in figure 3c to estimate the number of rare variants that are within 20bp of splicing junctions
Script: [RV_near_junction_window.sh](./data_processing/genotypes/annotation/RV_near_junction_window.sh)
```
bash RV_near_junction_window.sh > log_file_name.txt 2>&1
```


[ref1](https://doi.org/10.1016/j.ajhg.2018.05.002) Ganna, A., Satterstrom, F. K., Zekavat, S. M., Das, I., Kurki, M. I., Churchhouse, C., … Neale, B. M. (2018). Quantifying the Impact of Rare and Ultra-rare Coding Variation across the Phenotypic Spectrum. American Journal of Human Genetics, 102(6), 1204–1211.  https://github.com/andgan/ultra_rare_pheno_spectrum

[ref2](https://doi.org/10.1038/nature19057) Lek, M., Karczewski, K. J., Minikel, E. V, Samocha, K. E., Banks, E., Fennell, T., … Exome Aggregation Consortium. (2016). Analysis of protein-coding genetic variation in 60,706 humans. Nature, 536(7616), 285–91. 



# 2. RNA-seq data processing

The following steps will allow you to go from adapter trimming to quantififcation.
If you want to run 2.1 to 2.4 with one script, you can use [RD_analysis.sh](./data_processing/RNAseq/RD_analysis.sh)

```
bash RD_analysis.sh <batch_number> > log_file_name.txt 2>&1 &
```

For step by step analysis, follow this:
## 2.1 Adapters trimming
Input: Fastq files

Output: Trimmed Fastq files
		
Parameters:

* `fastq_dir`: Input directory that contains FASTQ files. The following code assumes that there are gzipped, paired R1 and R2 files for each sample, and it excludes "Undetermined" FASTQ files generated by `bcl2fastq`.
* `cutadapt_script`: [trim_adapters_150bpreads.sh](./data_processing/RNAseq/trim_adapters_150bpreads.sh) #include path to script

```cd ${FASTQ_DIR}
ls *merge_R1.fastq.gz | sed 's/_/\t/'| awk '{print $1}' |awk -v fastq_dir=$FASTQ_DIR 'BEGIN{OFS="\t"}{print fastq_dir"/"$1"_merge_R1.fastq.gz",fastq_dir"/"$1"_merge_R1.trimmed.fastq.gz", fastq_dir"/"$1"_merge_R2.fastq.gz",fastq_dir"/"$1"_merge_R2.trimmed.fastq.gz", fastq_dir"/"$1"_merge_short_reads.fastq.gz",fastq_dir"/"$1"_fastq_trimadapters.log"}' | \
parallel --jobs 15 --col-sep "\t" "${cutadapt_script} {1} {2} {3} {4} {5} {6}" 
```

## 2.2. Alignment

Alignment of trimmed reads using STAR

Input: Trimmed fastq files

Output: Genome and transcriptome bam files

Parameters:

* `index_dir`: Directory containing STAR INDEX 
* `bam_dir`: Directory for BAM file output
* `STAR_script`:[STAR_alignment_rare_disease.sh] (./data_processing/RNAseq/STAR_alignment_rare_disease.sh) #include path to script
* `Reference genome`: [hg19.fa](http://hgdownload.cse.ucsc.edu/downloads.html#human)
* `Gene annotation File`: [Gencode v19 GTF] (https://www.encodeproject.org/files/gencode.v19.annotation/)



```
ls *merge_R1.trimmed.fastq.gz | sed 's/_/\t/'| awk '{print $1}' |awk -v fastq_dir=$FASTQ_DIR -v bam_dir=$BAM_DIR 'BEGIN{OFS="\t"}{print fastq_dir"/"$1"_merge_R1.trimmed.fastq.gz", fastq_dir"/"$1"_merge_R2.trimmed.fastq.gz",bam_dir"/"$1}' |\
	parallel --jobs 5 --col-sep "\t" "${STAR_script} {1} {2} {3} $<index_dir>/STAR_INDEX_OVERHANG_150/ 150"
```

## 2.3 Filters
### 2.3.1 Filter bam files for reads mapping uniquely and remove PCR duplicates
Input: Transcriptome bam

Output: Filtered Transcriptome bam

Parameters:
`FILTER_script`: [Filters_bam_uniq_mq30_duplicates_rare_disease_transcriptome.sh](./data_processing/RNAseq/Filters_bam_uniq_mq30_duplicates_rare_disease_transcriptome.sh) #include path to script

```
cd $BAM_DIR
ls ${BAM_DIR}/*.Aligned.toTranscriptome.out.bam | parallel --jobs 10 --col-sep "\t" "${FILTER_script} {1}"
```

### 2.3.2 Filter junction files for splicing analysis.

Filter junction file for junctions detected with at least 10 reads uniquely spanning.

Input: [SJ.out.tab](http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STAR.posix/doc/STARmanual.pdf) files for each sample. This file is an output of STAR alignement

Output: Filtered junction files (.SJ.out_filtered_uniq_$threshold.tab)

Parameters:

`Junc_file`: Junction file generated by STAR 
`threshold`: Number of unique reads spanning the junction (Set to 10 in the analysis)
`prefix`: prefix used for output (i.e sample name)


```
cat $Junc_file  | awk -v T=$threshold '$7>=T ' > $prefix.SJ.out_filtered_uniq_$threshold.tab
```

## 2.4 Quantification

Use RSEM v.1.2.21 for quantification to obtain gene-level and transcript-level quantification

Parameters:
`RSEM_script`: [rsem_calculate_expression_raredisease.sh](./data_processing/RNAseq/rsem_calculate_expression_raredisease.sh)  #include path to script

```
ls *.Aligned.toTranscriptome.out_mapq30_sorted_dedup_byread.bam | sed 's/\./\t/'|awk '{print $1}' |awk -v bam_dir=$BAM_DIR  'BEGIN{OFS="\t"}{print bam_dir"/"$1".Aligned.toTranscriptome.out_mapq30_sorted_dedup_byread.bam",$1}' |   \
	parallel   --jobs 10 --col-sep "\t" "bash ${RSEM_script} {1} {2}"
```

## 2.5 Gene expression data processing 
### 2.5.1 Expression normalization

Input: gene level expression levels as output by RSEM

Output: Matrix of normalized counts

Script: [sva_pipeline.r](./data_processing/expression_norm/sva_pipeline.r)

Steps: 

* Read in raw count data output by RSEM
* Filter genes for minimum expression threshold (in the paper - TPM>0.5 in >50% of samples from each batch)
* Log transform
* Scale and center to generate expression Z-scores
* Run SVA
* Regress out SVs and SV+regression splines for SVs significantly associated with sample batch
* Run QC checks

### 2.5.2 Expression Outlier analysis

Input: Corrected counts

Output: Gene Sample Zscores file

Script: [exp_outlier_count.r](./expression_outlier_analysis/exp_outlier_count.r)
* Read in corrected gene expression count data (corrected_counts/gene/blood/[gene_count].txt) and sample metadata ([metadata_file].txt)
* Subset samples to include blood samples only
* Center and scale gene expression counts matrix to generate Z-scores
* Define Z-score thresholds for expression outlier calling
* For each Z-score threshold, find number of outlier for each sample
* Compute signififance of difference between outlier counts for under-expression and over-expression outliers (hypergeometric test)
* Write output

## 2.6 Splicing data processing 

This set of scripts is generating junction ratios and calculating splicing Z-scores for every sample

### 2.6.1 Prerequisites
To get started you will need
* [Metadata](./metadata.md) file
* Tissue to perform the analysis on (i.e. Blood)
* Junction files generated during STAR alignment (here they have been filtered for junctions with at least 10 reads uniquely spanning)


### 2.6.2 Generating junction ratios
Parameters:
* `metadatafile.tsv`: file containing [metadata](./metadata.md) information regarding the samples to analyze.
* `tissue`: What tissue (from the metadata file) to analyze (i.e Blood)
* `junc_dir`: Directory containing filtered jcuntions obtained after 2.3.2
* `output_dir`: Directory in which output will be written
* `freeze_sample_only`: TRUE/FALSE statement wether to include only samples for which in_freeze column set to yes in [metadata](./metadata.md) file.
* `DGN_samples`: TRUE/FALSE statement wether to include or not DGN samples in the analysis
* `PIVUS_samples`:  TRUE/FALSE statement wether to include or not PIVUS samples in the analysis

Script: [splicing_outlier_analysis.sh](./splicing_outlier_analysis/splicing_outlier_analysis.sh)
Relies on [make_list_junction_files.R] (./splicing_outlier_analysis/make_list_junction_files.R) and [splicing_outlier.py](./splicing_outlier_analysis/splicing_outlier.py)

Example command line:
```
bash splicing_outlier_analysis.sh \
	metadatafile.tsv \ # metadata file to use
	Blood  \ # tissue of analysis
	$junc_dir \ # folder containing filtered junction files
	$output_dir \ # output folder
	TRUE \ # wether or not to perform on samples from the freeze only
	FALSE \ # wether or not to include DGN samples
	TRUE > log_file.txt 2>&1 & # wether or not to include PIVUS cohort in the analysis
```
### 2.6.3 Generating Z-scores from the ratios
This step impute missing splicing ratios and get Z-scores for splicing data.

Input: <prefix>_junc_ratios_filtered.txt from 2.6.2
	
Script: [splicing_ratio_to_zscores.R](./splicing_outlier_analysis/splicing_ratio_to_zscores.R)

```
Rscript splicing_ratio_to_zscores.R
```

# 3. Phenotypic information
Download latest version of Human Phenotype Ontology (https://hpo.jax.org)

## 3.1 Get phenotype to gene link from Human Phenotype Ontology:
[phenotype_to_gene.txt](http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastSuccessfulBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt)

## 3.2 Parse HPO database
* Download latest version of [HPO obo] (https://hpo.jax.org/app/download/ontology)
* Parse database to get parent and child terms [parse_hpoterms.py](./HPO/parse_hpoterms.py0

# 4. Combine information to highlight candidate genes

## 4.1 Expression outlier filtering
### 4.1.1 Combine expression outlier and variant information
Script [get_rare_var_gene.sh](./expression_outlier_analysis/get_rare_var_gene.sh)

This script takes the outlier files from [2.5.2](###2.5.2-expression-outlier-analysis) and maps rare variants (MAF < 0.01) in genes or +/- 10kb around genes. This is done on a per-sample level. The output is a file containing gene z-score, position of rare variant, allele frequency, and cadd for each sample-gene pair. The filename is "outliers_rare_var_combined_10kb.txt".

### 4.1.2 Filter expression outliers using genetic and phenotype information
This step is filtering expression outlier data according to different criteria.
Candidate genes at each filtering step are printed in output



Parameters:
* `metadata file` : [Metadata](./metadata.md) file containing smaple information
* `exac file`: containing pLI score across genes
* `HPO gene to phenotype file`
* `HPO phenotype to gene file`
* `HPO parent - child - alternative ID`: obtained with [parse_hpoterms.py](./HPO/parse_hpoterms.py)
* `outliers_rare_var_combined_10kb.txt`: expression Z-score results with RV information obtained in 4.1.1

Script: [filter_expression_candidates.R](./expression_outlier_analysis/filter_expression_candidates.R)


Filters:
* `EXP_OUTLIER`: under-expression outlier (Zscore<=-2)
* `EXP_OUTLIER_PLI`: under-expression outliers in genes with pLI score >=0.9
* `EXP_OUTLIER_HPO`: under-expression outliers in genes linked to the phenotype (HPO terms)
* `EXP_OUTLIER_RV`: under-expression outliers in genes with a rare variant within 10kb
* `EXP_OUTLIER_RV_CADD`: under-expression outliers in genes with a deleterious rare variant within 10kb (CADD score >=10)
* `EXP_OUTLIER_RV_CADD_HPO`:  under-expression outliers in genes linked to the phenotype (HPO terms) with a deleterious rare variant within 10kb (CADD score >=10).


## 4.2 Splicing outlier filtering
### 4.2.1 Cross splicing Z-score results with RV information

* distance=20bp
Script: [splicing_outlier_genes_with_RV_window.sh](./splicing_outlier_analysis/splicing_outlier_genes_with_RV_window.sh]
```
# sort splicing outlier file
sort -k1,1 -k2,2n  splicing_outlier_file.txt>  splicing_outlier_file_sorted.txt

# combine splicing outlier file with variant data within a 20bp window of junction files
bash splicing_outlier_genes_with_RV_window.sh  splicing_outlier_file_sorted.txt > <log_file> 2>&1 &
```

### 4.2.2 Filter splicing outliers using genetic and phenotype information
This step is filtering splicing outlier data according to different criteria.
Candidate genes at each filtering step are printed in output


Script: [filter_expression_candidates.R](./splicing_outlier_analysis/filter_expression_candidates.R)


Filters:

* `SPLI_OUTLIER`: Splicing outlier gene (at least one junction with |Z-score|>=2)
* `SPLI_OUTLIER_PLI`: Splicing outlier in genes with pLI score >=0.9
* `SPLI_OUTLIER_HPO`: Splicing outlier in genes linked to the phenotype (HPO terms)
* `SPLI_OUTLIER_RV`: Splicing outlier in genes with a rare variant within 20bp of the outlier junction
* `SPLI_OUTLIER_RV_CADD`: Splicing outlier in genes with a deleterious rare variant within 20bp of the outlier junction(CADD score >=10)
* `SPLI_OUTLIER_RV_CADD_HPO`:  Splicing outlier in genes linked to the phenotype (HPO terms) with a deleterious rare variant within 20bp of the outlier junction(CADD score >=10)


# 5. Highlight candidate genes
This last step consists in scrolling through candidates obtained both using expression outlier and splicing outlier information 
We recommend to use the most stringent filters as they show the best results at highlighting the causal gene in our analysis.
This most stringent filter consist in filtering the expression/splicing outliers for genes:
* with deleterious rare variants nearby
* that are phenotypically relevant

Those results correspond to files with prefix `EXP_OUTLIER_RV_CADD_HPO` and `SPLI_OUTLIER_RV_CADD_HPO`

If no candidate is prioritized using those lists, the user can then look at other filter steps.










