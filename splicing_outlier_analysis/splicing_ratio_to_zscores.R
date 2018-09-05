#!/bin/R
# Laure Fresard

# transform splicing ratios to Z-scores after imputing missing values.
# output - zscore file for each tested junction


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


source('/users/lfresard/repos/rare_disease/scripts/manuscript_analyses/Figures_source.R')

setwd('/srv/scratch/restricted/rare_diseases/analysis/outlier_analysis/splicing/for_freeze')

#--- MAIN
# Get metadata
print("read in metadata")
metadata = "/srv/scratch/restricted/rare_diseases/data/metadata/2018_06_12_Rare_Disease_Metadata.tsv"
metadata = read_tsv(metadata) %>% filter(in_freeze=="yes") %>% filter(source_of_RNA=="Blood")

# Get metadata info for pivus and dgn
meta_pivus="/srv/scratch/restricted/rare_diseases/data/metadata/PIVUS_RNASequencingInfo.csv"
meta_pivus=read_csv(meta_pivus) %>% filter(Age=="70") %>% filter(RunOK=="Yes")
meta_pivus2="/srv/scratch/restricted/rare_diseases/data/metadata/Pivus_expCovariatesAll_unscaled.txt"
meta_pivus2=read_tsv(meta_pivus2)%>% filter(Age==0) #%>% filter(RunOK=="Yes")


# select samples for which we have genetic data
sample_with_data=metadata %>% filter(variant_data=="exome" | variant_data=="genome")%>% select(sample_id) %>% pull


# Get annotated junctions
print("read in annotated junctions")

annot_junc="/srv/scratch/restricted/rare_diseases/data/splicing/juncfiles/gencodev19_intronsstartplus1_proteincoding.genenames_uniqjunc.tsv"
annot_junc=read_tsv(annot_junc, col_names=FALSE)
colnames(annot_junc)=c("chr", "junc_start", "junc_end", "gene")
annot_junctions=paste(annot_junc$chr, annot_junc$junc_start, annot_junc$junc_end, annot_junc$gene, sep="_")



## affected status data
print("read in affected status")

affected_status_df=read.table('/srv/scratch/restricted/rare_diseases/analysis/outlier_analysis/splicing/for_freeze/ifrun9/sample_affected_status_freeze_RD_PIVUS.tsv', header=F)
colnames(affected_status_df)=c('sample', 'status')

affected_status_df$cohort=c(rep("RD",87), rep("PIVUS", 68))

samples=affected_status_df$sample
RD_samples=metadata$sample_id



## get ratio data
#ratio_file="/srv/scratch/restricted/rare_diseases/analysis/outlier_analysis/splicing/for_freeze/RD_freeze_junc_ratios.txt"

print("read in ratio data")

ratio_file="/srv/scratch/restricted/rare_diseases/analysis/outlier_analysis/splicing/for_freeze/ifrun9/RD_PIVUS_freeze_junc_ratios.txt"
ratios=as.data.frame(read.table(ratio_file,  header=T, sep="\t"))

all_samples=c(metadata$sample_id, meta_pivus$RNAseq_ID)


ratios=ratios %>% select(-Cluster) %>% distinct
# make junction names out of coordinates
junctions=paste(ratios$chr,ratios$junction_start, ratios$junction_end, ratios$gene, sep="_")
#samples_batch8=metadata%>% filter(batch==8)%>% select(sample_id) %>% pull


ratio.df=data.frame( ratios[,all_samples])
rownames(ratio.df)=junctions
ratio.df[is.na(ratio.df)]<-NA


print("remove junctions for which we have less than 30 samples with any value")

# remove junctiosn for which we have less than 30 samples with any value
ratio.df=ratio.df[colSums(!is.na(t(ratio.df)))>30,]

na.data=is.na(ratio.df)


print("impute missing data")

res.comp = imputePCA(t(ratio.df),ncp=10)
pca=prcomp(res.comp$completeObs)



sum_pca=t(summary(pca)$importance)
colnames(sum_pca)=c("sd", "prop_var", "cumul_var")
pcs_number=as.data.frame(sum_pca) %>% filter(cumul_var <=0.95) %>% nrow
pcs_selec=pca$x[,1:pcs_number]
mod = model.matrix(~1, data=as.data.frame(res.comp$completeObs))
#mod = model.matrix(~1, data=as.data.frame(t(ratio.df)))
modsv <- cbind(mod, pcs_selec)
fitsv <- lm.fit(modsv, res.comp$completeObs)


junction_zscore_matrix=as.matrix(fitsv$residuals)
junction_zscore_matrix.t=t(junction_zscore_matrix)
junction_zscore_matrix.t[na.data]<-NA
junction_zscore_matrix=t(junction_zscore_matrix.t)

print("scale and center data")

## Scale and center
junction_zscore_scale =scale(junction_zscore_matrix, center=TRUE, scale=TRUE)

print("format results")

junction_zscore_scale.t=t(junction_zscore_scale)
junction_zscore_scale.t.df=as.data.frame(junction_zscore_scale.t)
junction_zscore_scale.t.df$annotation_status=sapply(rownames(junction_zscore_scale.t.df),is_annotated, junc_ref_vect=annot_junctions)
junction_zscore_scale.t.df$chr=unlist(strsplit(rownames(junction_zscore_scale.t.df), "_"))[c(TRUE,FALSE,FALSE,FALSE)]
junction_zscore_scale.t.df$junc_start=unlist(strsplit(rownames(junction_zscore_scale.t.df), "_"))[c(FALSE,TRUE,FALSE,FALSE)]
junction_zscore_scale.t.df$junc_end=unlist(strsplit(rownames(junction_zscore_scale.t.df), "_"))[c(FALSE,FALSE,TRUE,FALSE)]
junction_zscore_scale.t.df$gene=unlist(strsplit(rownames(junction_zscore_scale.t.df), "_"))[c(FALSE,FALSE,FALSE, TRUE)]


junction_zscore_scale.t.df.m=melt(junction_zscore_scale.t.df, id.var=c("chr", "junc_start", "junc_end","gene", "annotation_status"))
junction_zscore_scale.t.df.m$absZ=abs(junction_zscore_scale.t.df.m$value)
zscores_blood=junction_zscore_scale.t.df.m %>% left_join(affected_status_df, by=c("variable"="sample")) %>% select(chr ,junc_start  ,junc_end, gene ,annotation_status ,variable, value, absZ, status)

print("write down results")

write.table(zscores_blood,'RD_freeze_junc_outliers_PCAratio_blood_nf_PIVUS_withPIVsamples_withmisdata.txt', col.names=T, row.names=F, quote=F, sep="\t", append=F)






