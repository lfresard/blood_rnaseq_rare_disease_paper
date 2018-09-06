#!/bin/R

# Run logistic regression model as shown in Fig 2A

## Load packages
library(data.table)
library(annotables)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(readr)
library(plyr)

setwd("[path/to/working/directory/]") # set working directory

##------------------- FUNCTIONS

## Function to get all outliers at a specified percentile threshold
pick_outliers <- function(medz, quantile_thresh, method) {
	if (method == "<") return(apply(t(medz), 1, function(x) {
			get_quantile <- as.numeric(quantile(x, quantile_thresh))
		 	list(which(x <= get_quantile))
		 }))
	else if (method == ">") return(apply(t(medz), 1, function(x) {
			get_quantile <- as.numeric(quantile(x, quantile_thresh))
		 	list(which(x >= get_quantile))
		 }))
	else print("Error: choose '<' or '>'")
}

## Function to return all outliers passing a given abs(Z-score) threshold
pick_global_outliers <- function(medz, thresh) {
	return(apply(t(medz), 1, function(x) list(which(abs(x) > thresh))))
}

## Function to count number of outliers for cases and controls
count_outlier_ids <- function(outlier_row, ids_filter, cases, controls) {
	outliers_in_gene <- ids_filter$sample_id[as.numeric(unlist(outlier_row))] # sample IDs of outliers
	sum_case <- length(which(outliers_in_gene %in% cases))
	sum_control <- length(which(outliers_in_gene %in% controls))
	return(c(sum_case, sum_control, sum_case/(sum_case + sum_control)))
}

## Function to return numeric p-value in asterick notation
emp_ast <- function(pvalues) {
	sapply(pvalues, function(i) {
		if (i > 0.05) return("")
		else if (i <= 0.5 & i > 0.01) return("*")
		else if (i <= 0.01 & i > 0.001) return("**")
		else return("***")
		})
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
stopifnot(ids$sample_id==gene_zscore[, 1]) # sample IDs must match exactly

## Drop rownames and make matrix numeric
gene_zscore_mat <- as.matrix(gene_zscore[, -1]) # convert dataframe to matrix
gene_zscore_num <- apply(gene_zscore_mat, 2, as.numeric) # convert to numeric
rownames(gene_zscore_num) <- unlist(gene_zscore[, 1]) # add back rownames
colnames(gene_zscore_num) <- unlist(strsplit(as.character(colnames(gene_zscore_num)), '[.]'))[c(TRUE,FALSE)] # Strip suffix from gene IDs in columns

## Scale and center
gene_zscore_scale <- scale(gene_zscore_num, center=TRUE, scale=TRUE)

## Save gene names to vector
genes <- colnames(gene_zscore_num)

## Get constraint metrics for current gene set
constraint_metrics_all <- read.table("[ExAC_constraint_metrics].txt", sep="\t", header=TRUE, stringsAsFactor=FALSE)
constraint_metrics <- merge(constraint_metrics_all[, c("gene", "syn_z", "mis_z", "lof_z", "pLI")], grch37[, c("ensgene", "symbol")], by.x="gene", by.y="symbol", all=FALSE)
constraint_metrics <- constraint_metrics[-which(duplicated(constraint_metrics)), ] # drop duplicated rows
constraint_metrics <- constraint_metrics[, c("ensgene", "pLI", "mis_z", "syn_z", "lof_z")] # re-arrange columns
constraint_metrics <- constraint_metrics[which(constraint_metrics$ensgene %in% genes), ] # filter to include only genes in current corrected matrix

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

## Get counts for case and control
case <- ids_filter$sample_id[ids_filter$affected_status == "Case"]
control <- ids_filter$sample_id[ids_filter$affected_status == "Control"]

## Set percentile and global outlier thresholds
iterator <- c(0.5, 0.25, 0.1, (10/nrow(ids_filter)), (5/nrow(ids_filter)), (2/nrow(ids_filter)), (1/nrow(ids_filter))) # percentile thresholds
quantile_method_array <- c(rep("<", length(iterator)), rep(">", length(iterator)-1)) # set which direction of distribution to look at at each percentile threshold
iterator <- unique(c(iterator, 1-iterator))
iterator <- iterator[order(iterator)]
stopifnot(length(iterator) == length(quantile_method_array))

## Initialize array to collect logistic model output
plot_collect <- array(NA, c(length(iterator), 9, 3)) # collect output from logistic model at each threshold
colnames(plot_collect) <- c("quantile", "pvalue", "coefficient", "error_low", "error_high", "ngenes", "nsamples", "predictor", "ncase")

## Loop through percentile thresholds fitting logistic model
for (k in 1:length(iterator)) {
	## Set current percentile threshold and method
	quantile_thresh <- iterator[k]
	quantile_method <- quantile_method_array[k] 

	## Get expression outliers
	outlier_indices <- pick_outliers(gene_zscore_scale_filter, quantile_thresh, quantile_method)

	## Get outlier counts for each gene
	outlier_counts <- lapply(outlier_indices, function(outlier_current) count_outlier_ids(outlier_current, ids_filter, case, control))
	outlier_counts_matrix <- matrix(unlist(outlier_counts), nrow=length(outlier_counts), ncol=3, byrow=T)
	outlier_counts_df <- data.frame(outlier_counts_matrix)
	outlier_counts_df$GeneID <- names(outlier_indices)
	colnames(outlier_counts_df) <- c("CaseOutliers", "ControlOutliers", "ProportionCase", "GeneID")

	## Combine with constraint metrics
	constraint_metrics_merge <- merge(outlier_counts_df, constraint_metrics, by.x="GeneID", by.y="ensgene", all=FALSE)
	constraint_metrics_merge <- constraint_metrics_merge[, -1] # drop gene name

	## Remove zero counts (no outliers)
	constraint_metrics_merge <- constraint_metrics_merge[apply(constraint_metrics_merge, 1, function(row) sum(row[1:2])>0), ]

	## Logistic regression model
	pred <- c("lof_z", "mis_z", "syn_z")
	for (l in 1:length(pred)) {
		frm <- as.formula(paste("cbind(CaseOutliers, ControlOutliers) ~", pred[l]))
		model <- glm(formula=frm, data=constraint_metrics_merge, family=binomial)
		plot_collect[k,1,l] <- k # quantile threshold
		mod_summary <- summary(model)
		plot_collect[k,2,l] <- mod_summary[[12]][2, 4] # pvalue
		plot_collect[k,3,l] <- mod_summary[[12]][2, 1] # coefficient estimate
		plot_collect[k,4,l] <- mod_summary[[12]][2, 1] - (1.96*mod_summary[[12]][2, 2]) # CI lower
		plot_collect[k,5,l] <- mod_summary[[12]][2, 1] + (1.96*mod_summary[[12]][2, 2]) # CI upper
		plot_collect[k,6,l] <- nrow(constraint_metrics_merge)
		plot_collect[k,7,l] <- length(unlist(outlier_indices[[1]])) # number of samples chosen in each gene at current percentile threshold
		plot_collect[k,8,l] <- pred[l] # current ExAC predictor
		plot_collect[k,9,l] <- length(case) # number of cases
	}
}

## Plot logistic coefficients
plot_collect_df <- data.frame(apply(plot_collect, 2, function(for_col) for_col))
plot_collect_df[, 1:7] <- apply(plot_collect_df[, 1:7], 2, function(for_col) as.numeric(for_col))
plot_collect_df$predictor <- revalue(plot_collect_df$predictor, c("lof_z"="Loss of Function", "mis_z"="Missense", "syn_z"="Synonymous"))

# Write output
write.table(global_out, file="global_exp_outliers.txt", col.names=F, row.names=F)
write.table(plot_collect_df, file="logistic_model_output.txt", col.names=T, row.names=F)

