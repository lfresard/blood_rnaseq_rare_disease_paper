#!/bin/R

# percentile_zscore_mapping.r

## Load packages
library(data.table)
library(dplyr)
library(readr)
library(plyr)

setwd("[path/to/working/directory/]") # set working directory

##------------------- FUNCTIONS

## Function to return all outliers passing a given abs(Z-score) threshold
pick_global_outliers <- function(medz, thresh) {
	return(apply(t(medz), 1, function(for_row) list(which(abs(for_row) > thresh))))
}

##------------------- MAIN

## Set helpful variables
metadata_file <- "[metadata_file].txt"
corrected_data_file <- "corrected_counts/gene/blood/[gene_counts].txt" # corrected data file
medz_ind_filt <- 100 # thresholds for defining global outlier
global_zscore <- 4 # z score for defining outlier for global outlier count

## Read in data 
gene_zscore <- fread(corrected_data_file, sep=",", header=TRUE, stringsAsFactor=FALSE)

## Read metadata
rd <- read_tsv(metadata_file) %>% filter(in_freeze=="yes") %>% filter(source_of_RNA=="Blood") %>% select(sample_id, affected_status)
ld <- data.frame(as.character(unlist(gene_zscore[grep("LD", as.character(unlist(gene_zscore[, 1]))), 1])), "Control")
colnames(ld) <- c("sample_id", "affected_status")
ids <- rbind.data.frame(rd, ld, stringsAsFactors=F)
ids <- subset(ids, sample_id %in% as.character(unlist(gene_zscore[, 1])))
ids <- ids[match(as.character(unlist(gene_zscore[, 1])), as.character(unlist(ids$sample_id))), ] # match IDs if necessary
stopifnot(ids$sample_id==gene_zscore[, 1]) # sample IDs must much exactly

## Drop rownames and make matrix numeric
gene_zscore_mat <- as.matrix(gene_zscore[, -1]) # convert dataframe to matrix
gene_zscore_num <- apply(gene_zscore_mat, 2, as.numeric) # convert to numeric
rownames(gene_zscore_num) <- unlist(gene_zscore[, 1]) # add back rownames
colnames(gene_zscore_num) <- unlist(strsplit(as.character(colnames(gene_zscore_num)), '[.]'))[c(TRUE,FALSE)] # Strip suffix from gene IDs in columns

## Scale and center
gene_zscore_scale <- scale(gene_zscore_num, center=TRUE, scale=TRUE)

## Remove global expression outliers
outlier_indices <- pick_global_outliers(gene_zscore_scale, global_zscore) # get under-expression outliers
outlier_count <- array(NA, c(nrow(ids), 2))
outlier_count[, 1] <- ids$sample_id
outlier_unlist <- ids$sample_id[as.numeric(unlist(outlier_indices))]
for (i in 1:nrow(outlier_count)) outlier_count[i, 2] <- sum(outlier_unlist == outlier_count[i, 1])
global_out <- outlier_count[which(as.numeric(outlier_count[, 2]) > medz_ind_filt), 1] # list of global outliers

## Remove global expression outliers from corrected data
gene_zscore_scale_filter <- gene_zscore_scale[-which(rownames(gene_zscore_scale) %in% global_out), ] # remove global outliers from matrix
ids_filter <- ids[which(ids$sample_id %in% rownames(gene_zscore_scale_filter)), ] # update IDs to remove IDs not in filtered matrix

## Set iterator
iterator <- c(0.5, 0.25, 0.1, (10/nrow(ids_filter)), (5/nrow(ids_filter)), (2/nrow(ids_filter)), (1/nrow(ids_filter))) # percentile thresholds
iterator <- unique(c(iterator, 1-iterator))
iterator <- iterator[order(iterator)]

## Get percentile<->Z-score mapping across all genes
percentile_collect <- data.frame(apply(gene_zscore_scale_filter, 2, function(x) quantile(x, iterator[1])))
colnames(percentile_collect) <- "value"
percentile_collect$percentile <- 1

for (i in 2:length(iterator)) {
	message(i)
	tmp <- data.frame(apply(gene_zscore_scale_filter, 2, function(x) quantile(x, iterator[i])))
	colnames(tmp) <- "value"
	tmp$percentile <- i
	percentile_collect <- rbind.data.frame(percentile_collect, tmp)
}

write.table(percentile_collect, file="percentile_zscore_map.txt", col.names=T, row.names=F)
