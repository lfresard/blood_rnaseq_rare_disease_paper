#!/bin/R

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

# Master directory
dir = Sys.getenv('RARE_DIS_DIR')


# load required functions
source('Figures_source.R')


##--- MAIN

# List junction files 
junction_dir=paste(dir,"/data/splicing/juncfiles/all_unfiltjunc/",sep="")
junc_suffix=".SJ.out_gene_info.tsv"

# Read rare disease samples metadata
metadata = paste(dir,"/data/metadata/2018_06_12_Rare_Disease_Metadata.tsv",sep="")
metadata = read_tsv(metadata) %>% filter(in_freeze=='yes') %>% filter(source_of_RNA== "Blood")
RD_samples=metadata$sample_id

# Get metadata info for pivus and dgn
meta_pivus2=paste(dir,"/data/metadata/Pivus_expCovariatesAll_unscaled.txt",sep="")
meta_pivus2=read_tsv(meta_pivus2)%>% filter(Age==0) 
pivus_samples=meta_pivus2$Sample 

dgn_junc=list.files(path=junction_dir, pattern=glob2rx("LD*.SJ.out_gene_info.tsv"))
dgn_samples=str_extract(dgn_junc, "LD[0-9]+")


# get all samples id
samples=c(RD_samples,pivus_samples,dgn_samples)

# select junc files corresponding to samples
junc_files=paste0(junction_dir,samples,junc_suffix)
samples=samples[(!file.size(junc_files) == 0)]
junc_files=junc_files[(!file.size(junc_files) == 0)]

#make a dataframe containing all junction counts for selected samples
junc_info=lapply(samples,process_junc_file,junction_dir,junc_suffix)
merged.data.frame = Reduce(function(...) merge(...,  by="junction",all=T), junc_info)
colnames(merged.data.frame)=c("junction",samples)


# Get annotated junctions
annot_junc=paste(dir,"/data/splicing/juncfiles/gencodev19_intronsstartplus1_proteincoding.genenames_uniqjunc.tsv",sep="")
annot_junc=read_tsv(annot_junc, col_names=FALSE)
colnames(annot_junc)=c("chr", "junc_start", "junc_end", "gene")
annot_junctions=paste(annot_junc$chr, annot_junc$junc_start, annot_junc$junc_end, annot_junc$gene, sep="_")

# get number of annotated junction per gene
annot_junc_num_gene=melt(table(str_extract(annot_junc$gene, "ENSG[0-9]+")))

# Filter merged data frame for junctions that are annotated
annot.junc=merged.data.frame %>% filter(junction %in% annot_junctions)

# filter out junction for which less than 20% of samples with at least 5 reads.
annot.junc[is.na(annot.junc)]<-0 # put NA junctions to 0
annot.junc.filt=annot.junc[rowSums(annot.junc[,2:ncol(annot.junc)]>=5) >(length(samples) *20/100),]

# get genes for which we observe this coverage
nb_junction_cov=melt(table(str_extract(annot.junc.filt$junction, "ENSG[0-9]+")))
junc_genes=nb_junction_cov$Var1
nb_junction_cov=nb_junction_cov %>% left_join(annot_junc_num_gene, by="Var1")
names(nb_junction_cov)=c("gene", "covered", "annotated")
nb_junction_cov=nb_junction_cov%>% mutate(percent_covered=covered/annotated)



# Read in disease gene lists
disease_genes_dir=paste(dir,"/data/candidate_gene_lists/",sep="")
gene_lists=list.files(path=disease_genes_dir, pattern='genelist.txt')

disease_lists=lapply(gene_lists, read_in_disease_lists, path_to_file=disease_genes_dir)
names(disease_lists)=unlist(strsplit(gene_lists, "_genelist.txt"))

# get ensembl names for all gene lists
disease_lists_ens=lapply(disease_lists,check_gene_names)



# get overlap between genes with coverage and disease lists
disease_junc=data.frame(OMIM=length( disease_lists_ens$OMIM[disease_lists_ens$OMIM %in% junc_genes])/length(disease_lists_ens$OMIM), Neurology=length( disease_lists_ens$Neurology[disease_lists_ens$Neurology %in% junc_genes])/length(disease_lists_ens$Neurology), Hematology=length( disease_lists_ens$Hematology[disease_lists_ens$Hematology %in% junc_genes])/length(disease_lists_ens$Hematology), Ophthalmology=length( disease_lists_ens$Ophtalmology[disease_lists_ens$Ophtalmology %in% junc_genes])/length(disease_lists_ens$Ophtalmology))
disease_junc.m=melt(disease_junc)
disease_junc.m$variable=factor(disease_junc.m$variable, levels=c("Hematology", "Ophthalmology", "Neurology", "OMIM"))

save.image(file = paste(dir,"/data/FigureS2.in.RData",sep=""))

