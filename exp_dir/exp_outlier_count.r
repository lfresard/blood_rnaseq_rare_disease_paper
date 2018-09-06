#!/bin/R

# Script to find number of case and control outliers at different Z-score thresholds as shown in Fig 2C

library(data.table)
library(plyr)
library(ggplot2)
library(annotables)
library(RColorBrewer)
library(stringr)
library(dplyr)
library(readr)

setwd("[path/to/working/directory/]") # set working directory

##------------------- FUNCTIONS

## Function to return all outliers passing a given abs(Z-score) threshold
pick_global_outliers <- function(medz, thresh) {
	return(apply(t(medz), 1, function(for_row) list(which(abs(for_row) > thresh))))
}

## Get all outliers satisfying some z score threshold (return z score)
pick_outliers_abs_z <- function(medz, zscore_thresh) return(apply(t(medz), 1, function(for_row) list(for_row[abs(for_row) >= zscore_thresh])))

## Get all outliers satisfying some z score threshold (return sample index)
pick_outliers_z_id <- function(medz, zscore_thresh, method) {
	if (method == "<") return(apply(t(medz), 1, function(for_row) list(which(for_row < -zscore_thresh))))
	else if (method == ">") return(apply(t(medz), 1, function(for_row) list(which(for_row > zscore_thresh))))
	else print("Error: choose '<' or '>'")
}

## Count case outliers
count_outliers <- function(ids_in, outliers, cohort, direction) {
	ids_cohort <- ids_in$sample_id[grep(cohort, ids_in$affected_status)]
	outlier_ids <- unlist(strsplit(names(unlist(outliers)), '[.]'))[c(FALSE,TRUE)]
	count_table_cohort <- data.frame(ids_cohort, 0)
	colnames(count_table_cohort) <- c("sample_id", "freq")
	outlier_count_by_sample <- data.frame(table(outlier_ids))
	outlier_count_by_sample <- outlier_count_by_sample[outlier_count_by_sample$outlier_ids %in% ids_cohort, ]
	if (nrow(outlier_count_by_sample) > 0) {
		count_table_cohort$Direction <- direction
		count_table_cohort$freq[match(outlier_count_by_sample$outlier_ids, count_table_cohort$sample_id)] <- outlier_count_by_sample$Freq
	}
	return(count_table_cohort)
}

## Count genes
count_genes <- function(outliers) {
	unlist_genes <- unlist(strsplit(names(unlist(outliers)), '[.]'))[c(TRUE,FALSE)]
	return(length(unique(unlist_genes)))
}

##------------------- MAIN

## Set variables
metadata_file <- "[metadata_file].txt"
corrected_data_file <- "corrected_counts/gene/blood/[gene_counts].txt"
medz_ind_filt <- 500 # threshold for defining global outlier 

## Read in data 
gene_zscore <- fread(corrected_data_file, sep=",", header=TRUE, stringsAsFactor=FALSE)

## Read metadata
rd <- read_tsv(metadata_file) %>% filter(in_freeze=="yes") %>% select(sample_id, affected_status)
ld <- data.frame(as.character(unlist(gene_zscore[grep("LD", as.character(unlist(gene_zscore[, 1]))), 1])), "Control")
colnames(ld) <- c("sample_id", "affected_status")
ids <- rbind.data.frame(rd, ld, stringsAsFactors=F)
ids <- subset(ids, sample_id %in% as.character(unlist(gene_zscore[, 1])))
stopifnot(ids$sample_id==gene_zscore[, 1]) # sample IDs must much exactly

## Drop rownames and make matrix numeric
gene_zscore_mat <- as.matrix(gene_zscore[, -1]) # convert dataframe to matrix
gene_zscore_num <- apply(gene_zscore_mat, 2, as.numeric) # convert to numeric
rownames(gene_zscore_num) <- unlist(gene_zscore[, 1]) # add back rownames
colnames(gene_zscore_num) <- unlist(strsplit(as.character(colnames(gene_zscore_num)), '[.]'))[c(TRUE,FALSE)] # Strip suffix from gene IDs in columns

## Scale and center
gene_zscore_scale <- scale(gene_zscore_num, center=TRUE, scale=TRUE)

## Save gene names to vector
genes <- colnames(gene_zscore_num)

## Remove global expression outliers
outlier_indices <- pick_global_outliers(gene_zscore_scale, 4) # get global expression outliers
outlier_count <- array(NA, c(nrow(ids), 2))
outlier_count[, 1] <- ids$sample_id
outlier_unlist <- ids$sample_id[as.numeric(unlist(outlier_indices))]
for (i in 1:nrow(outlier_count)) outlier_count[i, 2] <- sum(outlier_unlist==outlier_count[i, 1])
global_out <- outlier_count[which(as.numeric(outlier_count[, 2]) > medz_ind_filt), 1] # list of global outliers

## Remove global expression outliers from corrected data
gene_zscore_scale_filter <- gene_zscore_scale[-which(rownames(gene_zscore_scale) %in% global_out), ] # remove global outliers from matrix
ids_filter <- ids[which(ids$sample_id %in% rownames(gene_zscore_scale_filter)), ] # update IDs to remove IDs not in filtered matrix

thresholds <- c(2,3,4,5)
plot_collect_list <- list()
hyper_test_list <- list()

for (i in thresholds) {
	outlier_under <- pick_outliers_z_id(gene_zscore_scale_filter, i, "<") # get under-expression outliers
	outlier_over <- pick_outliers_z_id(gene_zscore_scale_filter, i, ">") # get over-expression outliers

	# unlist and combine
	for (k in 1:length(outlier_under)) {
		outlier_under[[k]] <- unlist(outlier_under[[k]])
		outlier_over[[k]] <- unlist(outlier_over[[k]])
	}

	# get outlier counts for case and control
	for (j in c("Case", "Control")) {
		num_case_low <- count_outliers(ids_filter, outlier_under, j, "Under")
		num_case_high <- count_outliers(ids_filter, outlier_over, j, "Over")

		num_case_combine <- rbind.data.frame(num_case_low, num_case_high)
		if (nrow(num_case_combine) > 0) {
			num_case_combine$Cutoff <- i # z score cutoff of current iteration
			num_case_combine$Cohort <- j # cohort
			num_case_combine$NGenes <- ifelse(num_case_combine[, 3]=="Under", 
							count_genes(outlier_under), count_genes(outlier_over)) # n genes at current threshold
		}
		plot_collect_list[[length(plot_collect_list)+1]] <- num_case_combine

		# Significance test (hypergeometric test)
		ids_success <- filter(ids_filter, affected_status==j & grepl("RD", sample_id))
		directions <- c("Under", "Cohort")
		for (l in 1:length(directions)) {
			if (directions[l]=="Under") outliers <- outlier_under else outliers <- outlier_over
			hyp_x <- sum(count_outliers(ids_success, outliers, j, "-")[, 2]) 
			hyp_m <- nrow(ids_success) * ncol(gene_zscore_scale_filter) 
			hyp_n <- (nrow(ids_filter) * ncol(gene_zscore_scale_filter)) - hyp_m 
			hyp_k <- length(unlist(outliers)) 
			hyp_prob <- phyper(hyp_x, hyp_m, hyp_n, hyp_k, lower.tail=F)
			hyper_test_write <- data.frame(ifelse(directions[l]=="Under", -i, i), hyp_prob, j)
			colnames(hyper_test_write) <- c("zscore", "pvalue", "cohort")
			hyper_test_list[[length(hyper_test_list)+1]] <- hyper_test_write
		}
	}
}

## Formatting
plot_collect <- rbind.fill(plot_collect_list)
colnames(plot_collect) <- c("SampleID", "OutlierCount", "Direction", "Cutoff", "Cohort", "NumGenes")
hyper_test_collect <- rbind.fill(hyper_test_list)
colnames(hyper_test_collect) <- c("zscore", "pvalue", "cohort")

## Reorder factor columns
plot_collect$Cutoff <- factor(plot_collect$Cutoff)
plot_collect$Direction <- factor(plot_collect$Direction, levels=c("Under", "Over"))
plot_collect$Direction <- revalue(plot_collect$Direction, c("Under"="Under Expression", "Over"="Over Expression"))

## Number of outliers per individual
p1 <- ggplot(plot_collect) + 
	geom_boxplot(aes(x=Cutoff, y=OutlierCount+1, group=interaction(Cutoff, Direction), colour=Direction),
	 width=0.8, position=position_dodge(width=0.85)) +
	labs(x="Z-score", y="Number of Outlier Genes") +
	theme_bw() +
	theme(strip.background=element_blank(),
		panel.grid.major=element_blank(), 
		panel.grid.minor=element_blank(),
		axis.text=element_text(size=10),
		panel.border=element_blank(),
		axis.line=element_line(size=0.5)) +
	scale_y_log10(breaks=c(0,1,10,100,500,2000)) +
	scale_colour_manual(values=c("gray40", "gray1"))
ggsave("outlier_count.pdf", p1, width=7, height=3.5)




