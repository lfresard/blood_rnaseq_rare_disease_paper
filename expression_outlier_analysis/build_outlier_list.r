## Rare variant enrichment analysis
## Code author: Craig Smail

library(annotables)
library(ggplot2)
library(plyr)
library(data.table)
library(parallel)
library(dplyr)
library(readr)
library(stringr)
library(ggsignif)

setwd("/srv/scratch/restricted/rare_diseases/analysis/outlier_analysis/expression_level/") 

args <- commandArgs(trailingOnly=TRUE)
corrected_data_file <- as.character(args[1])
metadata_file <- as.character(args[2])

## <------------------- FUNCTIONS

## Get all outliers satisfying some z score threshold (return z score)
pick_outliers_abs_z <- function(medz, zscore_thresh) {
	return(apply(t(medz), 1, function(for_row) list(for_row[abs(for_row) >= zscore_thresh])))
}

## <------------------- MAIN

## Read in data
metadata <- read_tsv(metadata_file)

## IDs in desired cohort
ids <- metadata %>% filter(source_of_RNA=="Blood") %>% filter(in_freeze=="yes") %>% filter(variant_data=="genome" | variant_data=="exome")

## Read and process corrected expression counts
gene_zscore <- fread(corrected_data_file, sep=",", header=TRUE, stringsAsFactor=FALSE)
gene_zscore_mat <- as.matrix(gene_zscore[, -1]) # convert dataframe to matrix
gene_zscore_num <- apply(gene_zscore_mat, 2, as.numeric) # convert to numeric
rownames(gene_zscore_num) <- unlist(gene_zscore[, 1]) # add back rownames
colnames(gene_zscore_num) <- unlist(strsplit(as.character(colnames(gene_zscore_num)), '[.]'))[c(TRUE,FALSE)] # Strip suffix from gene IDs in columns
gene_zscore_scale <- scale(gene_zscore_num, center=TRUE, scale=TRUE)
gene_zscore_scale <- gene_zscore_scale[which(rownames(gene_zscore_scale) %in% ids$sample_id), ]

## Get outliers
outliers <- pick_outliers_abs_z(gene_zscore_scale, 0)
outlier_list <- data.frame(unlist(outliers))
outlier_list <- cbind.data.frame(str_split_fixed(rownames(outlier_list), "\\.", 2), outlier_list, stringsAsFactors=F)
colnames(outlier_list) <- c("gene", "sample_id", "zscore")
outlier_list <- subset(outlier_list, sample_id %in% ids$sample_id)

## Merge with genomic information
outlier_list_merge <- merge(outlier_list, grch37, by.x="gene", by.y="ensgene", all=FALSE)
outlier_list_merge <- outlier_list_merge %>% select(chr, start, end, zscore, sample_id, gene)
chrom_order <- c((1:22),"X","Y","M")
outlier_list_merge$chr <- factor(outlier_list_merge$chr, levels=chrom_order, ordered=TRUE)
outlier_list_merge <- outlier_list_merge[do.call(order, outlier_list_merge[, c("chr", "start", "end", "zscore", "sample_id")]), ]

## Write to file by sample ID
for (i in unique(outlier_list_merge$sample_id)) {
	tmp <- subset(outlier_list_merge, sample_id==i)
	write.table(tmp, paste0("outliers/", i, "_outliers.bed"), row.names=F, col.names=F, sep="\t", quote=F)
}

