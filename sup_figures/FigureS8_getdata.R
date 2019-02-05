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



# plot the percent of sample(with gen data) for which there is at least one candidate accross filters

# Master directory
dir = Sys.getenv('RARE_DIS_DIR')



# Main
rm(list=ls())
load(file = paste0(dir,"/analysis/manuscript/figures_revision/Figure2.in.RData"))
load(file = paste0(dir,"/analysis/manuscript/figures_revision/Figure3.in.RData"))


# Get metadata
metadata = paste(dir,"/data/metadata/2018_12_02_Rare_Disease_Metadata.tsv",sep="")
metadata = read_tsv(metadata) %>% filter(in_freeze=="yes") %>% filter(source_of_RNA=="Blood", sequencing_status=="PASSED")
sample_with_data=metadata %>% filter(variant_data=="exome" | variant_data=="genome")%>% select(sample_id) %>% pull


# Expression non zero individual in %
#non_zero_samples=colSums(sample_exp_outlier_filter_withsgl_10kb_up.df[,5:(ncol(sample_exp_outlier_filter_withsgl_10kb_up.df))] != 0)/length(sample_with_data) *100
#non_zero_samples=melt(non_zero_samples)
#non_zero_samples$filter=c(1:13)
#

percent_expoutlier=data.frame(non_zero=exp_outlier_number.df %>% group_by(filter) %>% filter(value!=0) %>% summarize(percent_sample=length(sample)/length(cases_withGen)*100)%>% select(percent_sample)%>%pull,
	over5=exp_outlier_number.df %>% group_by(filter) %>% filter(value>5) %>% summarize(percent_sample=length(sample)/length(cases_withGen)*100)%>% select(percent_sample)%>%pull,
	over25=exp_outlier_number.df %>% group_by(filter) %>% filter(value>25) %>% summarize(percent_sample=length(sample)/length(cases_withGen)*100) %>% select(percent_sample)%>%pull, 
	#over50=exp_outlier_number.df %>% group_by(filter) %>% filter(value>50) %>% summarize(percent_sample=length(sample)/length(cases_withGen)*100)%>% select(percent_sample)%>%pull,
	#over100=exp_outlier_number.df %>% group_by(filter) %>% filter(value>100) %>% summarize(percent_sample=length(sample)/length(cases_withGen)*100)%>% select(percent_sample)%>%pull,
	filter=exp_outlier_number.df %>% group_by(filter) %>% filter(value>5) %>% summarize(percent_sample=length(sample)/length(cases_withGen)*100)%>% select(filter))

percent_expoutlier$filter=factor(percent_expoutlier$filter,levels=c("RV_10KB", "RV_10KB_CADD", "RV_10KB_HPO", "RV_10KB_CADD_HPO", "HPO", "EXP_OUTLIER", "EXP_OUTLIER_ASE", "EXP_OUTLIER_PLI", "EXP_OUTLIER_HPO","EXP_OUTLIER_RIVER", "EXP_OUTLIER_RV", "EXP_OUTLIER_RV_CADD", "EXP_OUTLIER_RV_CADD_HPO"))

percent_expoutlier.m=melt(percent_expoutlier, id.var=c("filter"))



#sample_exp_outlier_filter_withsgl_10kb_up.df.cases=sample_exp_outlier_filter_withsgl_10kb_up.df%>% filter(affected_status=="Case", sample_id %in% cases_withGen)
#non_zero_samples=colSums(sample_exp_outlier_filter_withsgl_10kb_up.df.cases[,5:(ncol(sample_exp_outlier_filter_withsgl_10kb_up.df.cases))] != 0)/length(cases_withGen) *100
#non_zero_samples=melt(non_zero_samples)
#non_zero_samples$filter=c(1:13)

# Splicing non zero individual in %
percent_splioutlier=data.frame(non_zero=res_splicing_number.df %>% group_by(variable) %>% filter(value!=0) %>% summarize(percent_sample=length(SAMPLE)/length(cases_withGen)*100)%>% select(percent_sample)%>%pull,
	over5= res_splicing_number.df %>% group_by(variable) %>% filter(value>5) %>% summarize(percent_sample=length(SAMPLE)/length(cases_withGen)*100)%>% select(percent_sample)%>%pull,
	#over25=res_splicing_number.df %>% group_by(variable) %>% filter(value>25) %>% summarize(percent_sample=length(SAMPLE)/length(cases_withGen)*100) %>% select(percent_sample)%>%pull, 
	#over50=res_splicing_number.df %>% group_by(filter) %>% filter(value>50) %>% summarize(percent_sample=length(sample)/length(cases_withGen)*100)%>% select(percent_sample)%>%pull,
	filter=res_splicing_number.df %>% group_by(variable) %>% filter(value>5) %>% summarize(percent_sample=length(SAMPLE)/length(cases_withGen)*100)%>% select(variable))
colnames(percent_splioutlier)=c("non_zero", "over5", "filter")
percent_splioutlier$filter=factor(percent_splioutlier$filter,levels=c("RV", "RV_CADD", "RV_HPO", "RV_CADD_HPO", "RV_JUNC", "RV_JUNC_CADD","RV_JUNC_HPO","RV_JUNC_CADD_HPO","HPO", "SPLI_OUTLIER", "SPLI_OUTLIER_PLI", "SPLI_OUTLIER_HPO", "SPLI_OUTLIER_RV","SPLI_OUTLIER_RV_CADD", "SPLI_OUTLIER_RV_CADD_HPO"))
percent_splioutlier.m=melt(percent_splioutlier, id.var=c("filter"))


#spli_non_zero_samples=colSums(sample_outlier_filter_withsgl.df[,3:(ncol(sample_outlier_filter_withsgl.df)-3)] != 0)/length(sample_with_data) *100
#spli_non_zero_samples=melt(spli_non_zero_samples)
#spli_non_zero_samples$filter=c(1:8)

save(metadata,sample_with_data,percent_expoutlier,percent_splioutlier,file = paste0(dir,"/analysis/manuscript/figures_revision/FigureS8.in.RData"))
