#!/bin/R

#Libraries
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(sva)
library(cowplot)
library(purrr)
library(broom)
library(reshape2)
library(ggfortify)
library(gridExtra)
library(data.table)
library(annotables)
library(corrplot)
library(RColorBrewer)
library(qvalue) 
library(ggpubr)
library(missMDA)


rm(list=ls())

# Master directory
dir = Sys.getenv('RARE_DIS_DIR')

# Load neessary functions
source(paste(dir,'/scripts/manuscript_analyses/Figures_source.R')


# Main
zscores_blood =as.data.frame(read.table(paste(dir,"/analysis/outlier_analysis/splicing/for_freeze/RD_freeze_junc_outliers_PCAratio_blood_nf_PIVUS_withPIVsamples_withmisdata.txt", sep=""), sep="\t", header=T))

zscores_blood =zscores_blood %>% mutate(gene=unlist(strsplit(as.character(gene),"[.]"))[c(TRUE,FALSE)]) 
genes=zscores_blood %>% select(gene) %>% distinct %>% pull



## Analyze results filtered for Rare Variants
RV_outlier_file=paste(dir,"/analysis/outlier_analysis/splicing/for_freeze/RV/RD_freeze_junc_outliers_PCAratio_blood_nf_PIVUS_withPIVsamples_withmisdata_sorted_RV_window.txt",sep="")
RV_outlier=as.data.frame(read.table(RV_outlier_file, sep="\t", header=F))
colnames(RV_outlier)=c(colnames(zscores_blood), "chr_var", "start_var", "end_var", "ref_var", "alt_var", "gnomAD_AF", "raw_cadd", "phred_cadd", "chr_gene", "start_gene", "end_gene", "gene_name")

# make a column with singleton included or expluded form the analysis
RV_outlier =RV_outlier %>% mutate(gnomAD_AF_withsgl = replace(gnomAD_AF, gnomAD_AF==".", 0))
RV_outlier =RV_outlier %>% mutate(gnomAD_AF_nosgl = replace(gnomAD_AF, gnomAD_AF==".", NA))


# format variant data so that only the max observed allele frequency is kept for variant were several frequencies are available
RV_outlier$gnomAD_AF_filt_withsgl=sapply(RV_outlier$gnomAD_AF_withsgl,process_variant)
RV_outlier$gnomAD_AF_filt_nosgl=sapply(RV_outlier$gnomAD_AF_nosgl,process_variant)

#filter for RV with MAF<=0.1%
maf=0.001
	# Excluding singletons
RV_outlier_filt_nosgl=RV_outlier %>% filter(!is.na(gnomAD_AF_filt_nosgl)) %>% filter(gnomAD_AF_filt_nosgl <=maf)
RV_outlier_filt_nosgl =RV_outlier_filt_nosgl %>% mutate(gene_name=str_extract(gene_name, "ENSG[0-9]+")) 
	# Including singletons
RV_outlier_filt_withsgl=RV_outlier %>% filter(!is.na(gnomAD_AF_filt_withsgl)) %>% filter(gnomAD_AF_filt_withsgl <=maf)
RV_outlier_filt_withsgl =RV_outlier_filt_withsgl %>% mutate(gene_name=str_extract(gene_name, "ENSG[0-9]+"))

RV_outlier_filt_withsgl_indels=RV_outlier_filt_withsgl %>% filter(nchar(as.character(alt_var))>1 | nchar(as.character(ref_var))>1)
RV_outlier_filt_withsgl_snps=RV_outlier_filt_withsgl %>% filter(nchar(as.character(alt_var))==1 & nchar(as.character(ref_var))==1)
# Make dataframe containing the number of genes with at least 1 outlier at different Zscore thresholds

sample_outlier.df=data.frame(sample_id=all_samples,
	"Z2"=sapply(all_samples, get_outlier_genes, threshold=2, outlier_df.m=zscores_blood, 0),
	"Z3"=sapply(all_samples, get_outlier_genes, threshold=3, outlier_df.m=zscores_blood, 0),
	"Z4"=sapply(all_samples, get_outlier_genes, threshold=4, outlier_df.m=zscores_blood, 0),
	"Z5"=sapply(all_samples, get_outlier_genes, threshold=5, outlier_df.m=zscores_blood, 0),
	"Z6"=sapply(all_samples, get_outlier_genes, threshold=6, outlier_df.m=zscores_blood, 0),
	affected_status=sapply(all_samples,get_affected_status, affected_status_df=affected_status_df))
	
# Add cohort information to data frame
sample_outlier.df=sample_outlier.df %>% left_join(affected_status_df, by=c("sample_id"="sample")) %>% select(-status) 

# transform data for plotting
sample_outlier.df.m=melt(sample_outlier.df,id.vars=c("sample_id","affected_status", "cohort"))
sample_outlier.df.m$variable=c(rep(2,152),rep(3,152),rep(4,152), rep(5,152), rep(6,152))

# Save data
save.image(file = paste(dir,"/data/FigureS11.in.RData",sep=""))
