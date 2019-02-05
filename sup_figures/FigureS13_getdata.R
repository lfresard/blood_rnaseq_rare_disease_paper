#!/bin/R

# Libraries
library(ggplot2)
library(dplyr)
library(reshape2)
library(cowplot)
library(readr)
library(RColorBrewer)
library(scales)


rm(list=ls())

# Master directory
dir = Sys.getenv('RARE_DIS_DIR')


# load required functions

source('Figures_source.R')
#source('/users/lfresard/repos/rare_disease/scripts/manuscript_analyses/Figures_source.R')



metadata = paste(dir,"/data/metadata/2018_12_02_Rare_Disease_Metadata.tsv",sep="")
metadata = read_tsv(metadata) %>% filter(in_freeze=='yes')  %>% filter(source_of_RNA=="Blood", sequencing_status=="PASSED")#%>% filter(variant_data=="exome" | variant_data=="genome" | variant_data=="none")


genome_data_summary=metadata %>% select(institution,affected_status,variant_data) %>% group_by(institution, affected_status,variant_data) %>%  summarize(variant_info_n=n())


# summary of varinat consequences per sample
variant_features="/users/xli6/projects/udn/feature_file/udn_variantcount.txt"
variant_features= read_tsv(variant_features)%>% filter(indiv_id %in% metadata$sample_id)
variant_features.lof= variant_features %>% select(indiv_id, variant_data,splice_count,frameshift_count,stop_count) %>% melt(id.vars=c("indiv_id", "variant_data") )


# Save data
save.image(file = paste(dir,"/analysis/manuscript/figures_revision/FigureS13.in.RData",sep=""))
