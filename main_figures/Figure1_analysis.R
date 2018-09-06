#!/bin/R
# LF
# November 2017


# This scripts is for figure 1 of the rare disease paper.
# Figure 1A: pie chart of broad disease categories
# Figure 1B: Expression in blood of disease genes
# Figure 1C: %genes expressed in blood when gene is expressed in only one tissue vs more than one tissue
# Figure 1D: expression of lof intolerant genes in blood

#--- Libraries
library(readr)
library(dplyr)
library(reshape2)
library(cowplot)
library(annotables)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(tximport)
library(stringr)
library(gtools)
library(ggpubr)


rm(list=ls())

source('/users/lfresard/repos/rare_disease/scripts/manuscript_analyses/Figures_source.R')

# Read rare disease samples metadata

metadata = "/srv/scratch/restricted/rare_diseases/data/metadata/2018_06_12_Rare_Disease_Metadata.tsv"
metadata = read_tsv(metadata) %>% filter(in_freeze=='yes') %>% filter(source_of_RNA== "Blood")

## Get number of case per disease category
categ_case=metadata %>% 
	filter(affected_status=="Case") %>% 
	group_by(disease_category) %>% 
	summarize(counts=n()) %>% 
	arrange(by=counts) 
categ_case$disease_category=c("Infectious diseases","Multiple Congenital \nAnomalies", "Rheumatology", "Allergies" ,"Cardiology", "Gynecology", "Other", "Musculoskeletal and \northopedics", "Opthalmology", "Hematology", "Neurology")

categ_case=categ_case %>% mutate(disease_category=factor(disease_category,levels=disease_category[length(disease_category):1]))

##-------------------------
## Expression of genes in blood by disease category
# Load Expression Data
gene_counts_RD = paste0("/srv/scratch/restricted/rare_diseases/data/quantification/rsem/",metadata$sample_id,".genes.results")
gene_counts_RD=gene_counts_RD[file.exists(gene_counts_RD)]
list_DGN_counts=list.files(path="/srv/scratch/restricted/rare_diseases/data/quantification/rsem/", pattern=glob2rx("LD*.genes.results*"))
genes_counts_DGN=paste0("/srv/scratch/restricted/rare_diseases/data/quantification/rsem/",list_DGN_counts)
DGN_names=str_extract(list_DGN_counts, "LD[0-9]+")
list_PIVUS_counts=list.files(path="/srv/scratch/restricted/rare_diseases/data/quantification/rsem/", pattern=glob2rx("PIVUS*.genes.results*"))[c(TRUE,FALSE)]
genes_counts_PIVUS=paste0("/srv/scratch/restricted/rare_diseases/data/quantification/rsem/",list_PIVUS_counts)
PIVUS_names=str_extract(list_PIVUS_counts, "PIVUS[0-9]+")

gene_counts=c(gene_counts_RD,genes_counts_DGN,genes_counts_PIVUS)


geneIdCol <- "gene_id"
abundanceCol <- "TPM"
countsCol <- "expected_count"
lengthCol <- "effective_length"
col.types <- readr::cols(
  readr::col_character(),readr::col_character(),readr::col_double(),readr::col_double(),
  readr::col_double(),readr::col_double(),readr::col_double()
)
importer <- function(x) readr::read_tsv(x, progress=FALSE, col_types=col.types)



# build TPM matrix
#rsem = tximport(gene_counts, type='rsem', txIn = FALSE, abundanceCol="TPM")
rsem = tximport(gene_counts, type='none', txIn = FALSE, abundanceCol=abundanceCol,geneIdCol=geneIdCol, countsCol=countsCol, lengthCol=lengthCol , importer=importer )

colnames(rsem$abundance) = c(metadata$sample_id,DGN_names,PIVUS_names)
rownames(rsem$abundance) = str_extract(rownames(rsem$abundance), "ENSG[0-9]+")
rsem$abundance = rsem$abundance[!is.na(rownames(rsem$abundance)),]


gene_ids = data.frame(ensgene = rownames(rsem$abundance))
tpms_df = tbl_df(cbind(gene_ids, rsem$abundance))


# Filter data for protein coding genes
protein_coding_genes=grch37 %>% filter(biotype=="protein_coding") %>% select(ensgene)
tpms_df_pc=tpms_df[ tpms_df$ensgene %in% protein_coding_genes$ensgene, ]


# Get percent of individuals with TPM < 1 for every protein coding gene
pct_tpm_below1=apply(tpms_df_pc[,2:ncol(tpms_df_pc)], 1, get_percent, threshold=1)


# Get different threshold of % individual expressed for every genes
pct_df=data.frame(GENE=tpms_df_pc$ensgene, 
	PCT_0= ifelse(pct_tpm_below1 <=100,1,0),
	PCT_10= ifelse(pct_tpm_below1 <=90,1,0),	
	PCT_20= ifelse(pct_tpm_below1 <=80,1,0),
	PCT_30= ifelse(pct_tpm_below1 <=70,1,0),
	PCT_40= ifelse(pct_tpm_below1 <=60,1,0),
	PCT_50= ifelse(pct_tpm_below1 <=50,1,0), 
	PCT_60= ifelse(pct_tpm_below1 <=40,1,0),
	PCT_70= ifelse(pct_tpm_below1 <=30,1,0), 
	PCT_80= ifelse(pct_tpm_below1 <=20,1,0),
	PCT_90= ifelse(pct_tpm_below1 <=10,1,0),
	PCT_100=ifelse(pct_tpm_below1 ==0,1,0))

# Read in disease gene lists
disease_genes_dir='/srv/scratch/restricted/rare_diseases/data/candidate_gene_lists/'
gene_lists=list.files(path=disease_genes_dir, pattern='genelist.txt')

disease_lists=lapply(gene_lists, read_in_disease_lists, path_to_file=disease_genes_dir)
names(disease_lists)=unlist(strsplit(gene_lists, "_genelist.txt"))

# get ensembl names for all gene lists
disease_lists_ens=lapply(disease_lists,check_gene_names)

# For each genes list, look at how many genes are expressed regarding to filters set above

pct_in_blood=lapply(disease_lists_ens,get_pct_disease_blood,pct_df=pct_df)
pct_in_blood_df=as.data.frame(do.call(rbind, pct_in_blood))
#pct_in_blood_df=pct_in_blood_df %>% select(PCT_20,PCT_50,PCT_90,DISEASE)
pct_in_blood_df$DISEASE=factor(rownames(pct_in_blood_df))
pct_in_blood_df =pct_in_blood_df %>% filter(DISEASE %in% c("Neurology", "Ophtalmology",  "OMIM", "Hematology"))
pct_in_blood_df$DISEASE=factor(pct_in_blood_df$DISEASE, levels=c("Neurology", "Hematology","Ophtalmology",  "OMIM"))
pct_in_blood_df.m=melt(pct_in_blood_df, by='DISEASE')
pct_in_blood_df.m$percent_ind=c(rep(0,4),rep(10,4), rep(20,4),rep(30,4),rep(40,4),rep(50,4),rep(60,4),rep(70,4),rep(80,4),rep(90,4),rep(100,4))



# Look at % genes expressed in blood in function of average TPM expression 

avg_tpm_pc=data.frame(GENE=tpms_df_pc$ensgene, avg_tpm=rowMeans(tpms_df_pc[,2:ncol(tpms_df_pc)]), median_tpm=apply(tpms_df_pc[,2:ncol(tpms_df_pc)], 1, median))


#breaks=c(0, 0.1,1,10,100, 100000)
breaks=c(0, 0.1,10,100000)
breaks2=c(0, 0.1,100000)
avg_tpm_pc=make_bins(avg_tpm_pc, breaks)
avg_tpm_pc2=make_bins(avg_tpm_pc, breaks2)

pct_in_blood_bin=lapply(disease_lists_ens,get_pct_disease_blood_bin,pct_df=avg_tpm_pc)
pct_in_blood_bin_df=as.data.frame(do.call(rbind, pct_in_blood_bin))
pct_in_blood_bin2=lapply(disease_lists_ens,get_pct_disease_blood_bin,pct_df=avg_tpm_pc2)
pct_in_blood_bin_df2=as.data.frame(do.call(rbind, pct_in_blood_bin2))

pct_in_blood_bin_df$DISEASE=unlist(strsplit(rownames(pct_in_blood_bin_df), "[.]"))[c(TRUE,FALSE)]
pct_in_blood_bin_df =pct_in_blood_bin_df %>% filter(DISEASE %in% c("Neurology", "Hematology","Ophtalmology",  "OMIM"))
pct_in_blood_bin_df2$DISEASE=unlist(strsplit(rownames(pct_in_blood_bin_df2), "[.]"))[c(TRUE,FALSE)]
pct_in_blood_bin_df2 =pct_in_blood_bin_df2 %>% filter(DISEASE %in% c("Neurology", "Hematology","Ophtalmology",  "OMIM"))


pct_in_blood_bin_df$DISEASE=factor(pct_in_blood_bin_df$DISEASE, levels=c("Neurology", "Hematology","Ophtalmology" , "OMIM"))
pct_in_blood_bin_df2$DISEASE=factor(pct_in_blood_bin_df2$DISEASE, levels=c("Neurology", "Hematology","Ophtalmology" , "OMIM"))


bin_df=data.frame(pct_in_blood_bin_df[order(pct_in_blood_bin_df$bin),] %>% select(bin) %>% distinct, INDEX=c(1:3))
bin_df2=data.frame(pct_in_blood_bin_df2[order(pct_in_blood_bin_df2$bin),] %>% select(bin) %>% distinct, INDEX=c(1,3))


pct_in_blood_bin_df$INDEX=sapply(pct_in_blood_bin_df$bin, get_index, bin_df=bin_df)
pct_in_blood_bin_df2$INDEX=sapply(pct_in_blood_bin_df2$bin, get_index, bin_df=bin_df2)
pct_in_blood_bin_df2=pct_in_blood_bin_df2 %>% filter(bin !="[0,0.1)")
#pct_in_blood_bin_df.m=melt(pct_in_blood_bin_df, by='DISEASE')

pct_in_blood_bin_df3=rbind(pct_in_blood_bin_df,pct_in_blood_bin_df2)%>% filter(bin !="[10,1e+05)")



# make equivalent analysis using median instead of average
median_tpm_pc=make_bins_median(avg_tpm_pc, breaks)
median_tpm_pc2=make_bins_median(avg_tpm_pc, breaks2)
pct_in_blood_bin_median=lapply(disease_lists_ens,get_pct_disease_blood_bin,pct_df=median_tpm_pc)
pct_in_blood_bin_median_df=as.data.frame(do.call(rbind, pct_in_blood_bin_median))
pct_in_blood_bin_median_df$DISEASE=unlist(strsplit(rownames(pct_in_blood_bin_median_df), "[.]"))[c(TRUE,FALSE)]
pct_in_blood_bin_median_df =pct_in_blood_bin_median_df %>% filter(DISEASE %in% c("Neurology", "Hematology","Ophtalmology",  "OMIM"))

pct_in_blood_bin_median2=lapply(disease_lists_ens,get_pct_disease_blood_bin,pct_df=median_tpm_pc2)
pct_in_blood_bin_median_df2=as.data.frame(do.call(rbind, pct_in_blood_bin_median2))
pct_in_blood_bin_median_df2$DISEASE=unlist(strsplit(rownames(pct_in_blood_bin_median_df2), "[.]"))[c(TRUE,FALSE)]
pct_in_blood_bin_median_df2 =pct_in_blood_bin_median_df2 %>% filter(DISEASE %in% c("Neurology", "Hematology","Ophtalmology",  "OMIM"))

bin_median_df=data.frame(pct_in_blood_bin_median_df[order(pct_in_blood_bin_median_df$bin),] %>% select(bin) %>% distinct, INDEX=c(1:3))
bin_median_df2=data.frame(pct_in_blood_bin_df2[order(pct_in_blood_bin_df2$bin),] %>% select(bin) %>% distinct, INDEX=c(1,3))



pct_in_blood_bin_median_df$INDEX=sapply(pct_in_blood_bin_median_df$bin, get_index, bin_df=bin_df)
pct_in_blood_bin_median_df2$INDEX=sapply(pct_in_blood_bin_median_df2$bin, get_index, bin_df=bin_df2)
pct_in_blood_bin_median_df2=pct_in_blood_bin_median_df2 %>% filter(bin !="[0,0.1)")

pct_in_blood_bin_median_df3=rbind(pct_in_blood_bin_median_df,pct_in_blood_bin_median_df2)%>% filter(bin !="[10,1e+05)")

pct_in_blood_bin_median_df3$DISEASE=factor(pct_in_blood_bin_median_df3$DISEASE, levels=c("OMIM","Neurology", "Ophtalmology", "Hematology"))

##-------------------------
# Genes expressed in more than 1 tissue and expression in Blood
# We looked at GTEx v7 RPKM for this plot

# Read in % results
exp_multtissues_file="/srv/scratch/restricted/rare_diseases/data/gtex_tissues_exp/Blood_genes_comp_othertissues_exp.txt"
exp_multtissues=read.table(exp_multtissues_file, header=T, sep="\t")


# Make data frame out of results
multis_df=data.frame(tissues=c('1 tissue', '>1 tissue'), 
	overlap_blood=c(exp_multtissues$gene_in_blood[exp_multtissues$tissue_number==1]*100, sum(exp_multtissues$gene_in_blood[2:nrow(exp_multtissues)])*100),
	N_genes=c(exp_multtissues$gene_number[exp_multtissues$tissue_number==1],sum(exp_multtissues$gene_number[2:nrow(exp_multtissues)])))

multis_df$tissues=factor(multis_df$tissues, levels=rev(levels(multis_df$tissues)))


## Look at gene level
exp_multtissues_genes_file="/srv/scratch/restricted/rare_diseases/data/gtex_tissues_exp/Number_of_tissue_expression.txt"
exp_multtissues_genes=as.data.frame(read.table(exp_multtissues_genes_file, header=F, sep="\t"))
colnames(exp_multtissues_genes)=c("ensgene", "tissues")
exac = as.data.frame(read.table('/users/lfresard/NumberOfIndividuals_Impact/enrichment_disease-conservation_database/annotations/EXAC/forweb_cleaned_exac_r03_march16_z_data_pLI.txt', sep = '\t', stringsAsFactors = F, header = T))

# stratify analysis between high/low pLIs
exp_multtissues_genes=exp_multtissues_genes %>% 
	inner_join(grch37, by="ensgene") %>% 
	select(ensgene, tissues, symbol) %>% 
	distinct %>% inner_join(exac, by=c("symbol"="gene")) %>% 
	select(ensgene,tissues,syn_z,mis_z,lof_z,pLI) %>% 
	distinct %>% mutate(multi=ifelse(tissues==1,"no", "yes"))%>%select(-one_of(c("tissues", "pLI")))


exp_multtissues_genes.m=melt(exp_multtissues_genes, id.vars = c("ensgene", "multi"))


## expression of extreme lofz and misz in blood

# get average expression in blood
tpms_df_avg=data.frame(ensgene=tpms_df$ensgene, avg_exp=rowMeans(tpms_df[,2:ncol(tpms_df)]))
tpms_df_avg=tpms_df_avg %>%
	inner_join(grch37, by="ensgene") %>% filter(biotype == "protein_coding")%>%
	select(ensgene, symbol,avg_exp) %>% 
	inner_join(exac, by=c("symbol"="gene")) %>% 
	select(ensgene,avg_exp,syn_z,mis_z,lof_z, pLI) %>% distinct %>% filter(pLI>=.9)

#breaks=c(0, 0.1,1,10,100,10000)
breaks=c(0, 0.1,10,10000,100000)
tpms_df_avg=transform(tpms_df_avg, bin = cut(avg_exp, breaks, right=F))




save.image(file="Figure1.RData")

