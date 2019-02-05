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

# Get metadata
metadata = "/srv/scratch/restricted/rare_diseases/data/metadata/2018_12_02_Rare_Disease_Metadata.tsv"
metadata = read_tsv(metadata) %>% filter(in_freeze=="yes",status=="PASSED", is_RD=="yes",source_of_RNA=="Blood")
# Get metadata info for pivus and dgn
meta_pivus="/srv/scratch/restricted/rare_diseases/data/metadata/PIVUS_RNASequencingInfo.csv"
meta_pivus=read_csv(meta_pivus) %>% filter(Age=="70") %>% filter(RunOK=="Yes")
meta_pivus2="/srv/scratch/restricted/rare_diseases/data/metadata/Pivus_expCovariatesAll_unscaled.txt"
meta_pivus2=read_tsv(meta_pivus2)%>% filter(Age==0) #%>% filter(RunOK=="Yes")


# get affected status data
affected_status_PIVUS=data.frame(sample_id=meta_pivus2$Sample, affected_status=rep("Control", nrow(meta_pivus2)))

affected_status_df=rbind(metadata%>% select(sample_id,affected_status),affected_status_PIVUS)
colnames(affected_status_df)=c('sample', 'status')

affected_status_df$cohort=c(rep("RD",nrow(metadata)), rep("PIVUS",  nrow(meta_pivus2)))
samples=affected_status_df$sample

# read in exac data
exac = as.data.frame(read.table('/users/lfresard/NumberOfIndividuals_Impact/enrichment_disease-conservation_database/annotations/EXAC/forweb_cleaned_exac_r03_march16_z_data_pLI.txt', sep = '\t', stringsAsFactors = F, header = T))
exac=exac %>% inner_join(grch37, by=c("gene"="symbol")) %>% select(ensgene,syn_z, mis_z, lof_z, pLI) 



zscores_blood =as.data.frame(read.table(paste(dir, "/analysis/outlier_analysis/splicing/2018_12_3_paper_revisions/2018_12_03_splicing_zscores.txt", sep=""), sep="\t", header=T))

zscores_blood =zscores_blood %>% mutate(gene=unlist(strsplit(as.character(gene),"[.]"))[c(TRUE,FALSE)]) 
genes=zscores_blood %>% select(gene) %>% distinct %>% pull

all_samples=c(metadata$sample_id, meta_pivus$RNAseq_ID)


## Analyze results filtered for Rare Variants
RV_outlier_file=paste(dir,"/analysis/outlier_analysis/splicing/2018_12_3_paper_revisions/RV/2018_12_03_splicing_zscores_sorted_RV_window.txt",sep="")
RV_outlier=as.data.frame(read.table(RV_outlier_file, sep="\t", header=F))
colnames(RV_outlier)=c(colnames(zscores_blood), "chr_var", "start_var", "end_var", "ref_var", "alt_var", "gnomAD_AF", "raw_cadd", "phred_cadd", "chr_gene", "start_gene", "end_gene", "gene_name")

# make a column with singleton included or expluded form the analysis
RV_outlier =RV_outlier %>% mutate(gnomAD_AF_withsgl = replace(gnomAD_AF, gnomAD_AF==".", 0))


# format variant data so that only the max observed allele frequency is kept for variant were several frequencies are available
RV_outlier$gnomAD_AF_filt_withsgl=sapply(RV_outlier$gnomAD_AF_withsgl,process_variant)

#filter for RV with MAF<=0.1%
maf=0.001
	# Including singletons
RV_outlier_filt_withsgl=RV_outlier %>% filter(!is.na(gnomAD_AF_filt_withsgl)) %>% filter(gnomAD_AF_filt_withsgl <=maf)
RV_outlier_filt_withsgl =RV_outlier_filt_withsgl %>% mutate(gene_name=str_extract(gene_name, "ENSG[0-9]+"))
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
sample_outlier.df.m$variable=c(rep(2,length(samples)),rep(3,length(samples)),rep(4,length(samples)), rep(5,length(samples)), rep(6,length(samples)))

# Save data
save.image(file = paste(dir,"/analysis/manuscript/figures_revision/FigureS11.in.RData",sep=""))
