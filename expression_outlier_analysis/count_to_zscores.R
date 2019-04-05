#!/bin/R

##Libraries
library(data.table)

## Read in data 
gene_zscore <- fread(corrected_data_file, sep=",", header=TRUE, stringsAsFactor=FALSE)


## Drop rownames and make matrix numeric
gene_zscore_mat <- as.matrix(gene_zscore[, -1]) # convert dataframe to matrix
gene_zscore_num <- apply(gene_zscore_mat, 2, as.numeric) # convert to numeric
rownames(gene_zscore_num) <- unlist(gene_zscore[, 1]) # add back rownames
colnames(gene_zscore_num) <- unlist(strsplit(as.character(colnames(gene_zscore_num)), '[.]'))[c(TRUE,FALSE)] # Strip suffix from gene IDs in columns



## Scale and center
gene_zscore_scale <- scale(gene_zscore_num, center=TRUE, scale=TRUE)

## Melt zscore matrix
gene_zscore_scale.m=melt(gene_zscore_scale)
colnames(gene_zscore_scale.m)=c("sample_id", "gene","zscore")

## Write down results

write.table(gene_zscore_scale.m, "outliers_zscore_pair_spline.txt", sep="\t", row.names=F, col.names=T, quote=F) 
