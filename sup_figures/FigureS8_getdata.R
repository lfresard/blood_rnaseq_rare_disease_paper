#!/bin/R

# plot the percent of sample(with gen data) for which there is at least one candidate accross filters

# Master directory
dir = Sys.getenv('RARE_DIS_DIR')


#Libraries
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(dplyr)
library(readr)
library(reshape2)

# Main
rm(list=ls())
load(file = paste(dir,"/analysis/outlier_analysis/expression_level/expression_filters_genes.RData"),sep="")
load(file = paste(dir,"/analysis/outlier_analysis/splicing/for_freeze/splicing_filters_genes.RData"), sep="")


# Get metadata
metadata = paste(dir,"/data/metadata/2018_06_12_Rare_Disease_Metadata.tsv",sep="")
metadata = read_tsv(metadata) %>% filter(in_freeze=="yes") %>% filter(source_of_RNA=="Blood")
sample_with_data=metadata %>% filter(variant_data=="exome" | variant_data=="genome")%>% select(sample_id) %>% pull


# Expression non zero individual in %
non_zero_samples=colSums(sample_exp_outlier_filter_withsgl_10kb_up.df[,3:(ncol(sample_exp_outlier_filter_withsgl_10kb_up.df)-3)] != 0)/length(sample_with_data) *100
non_zero_samples=melt(non_zero_samples)
non_zero_samples$filter=c(1:13)
 

# Splicing non zero individual in %
spli_non_zero_samples=colSums(sample_outlier_filter_withsgl.df[,3:(ncol(sample_outlier_filter_withsgl.df)-3)] != 0)/length(sample_with_data) *100
spli_non_zero_samples=melt(spli_non_zero_samples)
spli_non_zero_samples$filter=c(1:8)

save.image(file = paste(dir,"/data/FigureS8.in.RData",sep=""))
