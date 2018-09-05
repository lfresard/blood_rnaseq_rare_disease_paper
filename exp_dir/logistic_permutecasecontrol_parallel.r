#!/bin/R

# Run case/control label permutation script as shown in Sup Fig 2D

library(data.table)
library(annotables)
library(ggplot2)
library(gridExtra)
library(parallel)
library(readr)
library(dplyr)

## Function to return all outliers passing a given abs(Z-score) threshold
pick_global_outliers <- function(medz, thresh) {
	return(apply(t(medz), 1, function(for_row) list(which(abs(for_row) > thresh))))
}

## Function to get all outliers at a specified percentile threshold
pick_outliers <- function(medz, quantile_thresh) {
	if (quantile_thresh<=0.5) return(apply(t(medz), 1, function(for_row) {
			get_quantile <- as.numeric(quantile(for_row, quantile_thresh))
		 	list(which(for_row <= get_quantile))
		 }))
	else if (quantile_thresh>0.5 & quantile_thresh<=1) return(apply(t(medz), 1, function(for_row) {
			get_quantile <- as.numeric(quantile(for_row, quantile_thresh))
		 	list(which(for_row >= get_quantile))
		 }))
	else print("Percentile threshold must be between 0 and 1")
}

## Function to get labels for plotting
get_plot_labels <- function(n) {
	iterator <- c(0.5, 0.25, 0.1, 10/n, 5/n, 2/n, 1/n)
	nsamples <- round(n*iterator)
	nsamples <- c(nsamples[order(nsamples)], nsamples[-1])
	iterator <- unique(c(iterator, 1-iterator))
	iterator <- iterator[order(iterator)]
	for_return <- matrix(nrow=2, ncol=length(iterator))
	for_return[1, ] <- iterator
	for_return[2, ] <- nsamples
	return(for_return)
}

## Function to return a plot of results
ggplot_build <- function(dat, nsamples, iterator, emp) {
	p1 <- ggplot(dat) + 
	geom_hline(yintercept=0) + 
	geom_errorbar(aes(x=threshold, ymin=error_low, ymax=error_high), colour="gray") + 
	geom_point(aes(x=threshold, y=real_coefficient, fill=predictor), shape=21, colour="Black") +
	facet_grid(. ~ predictor, scales="free_y") +
	labs(x="Percentile", y="Log Odds") +
	scale_x_continuous(breaks=seq(1, length(nsamples), 1), labels=paste0(round(iterator*100, 3), "%")) +
	geom_text(aes(x=threshold, y=real_coefficient+0.01), label=emp, size=3) +
	RD_theme +
	theme(axis.text.x=element_text(angle=45, hjust=1),
		strip.text=element_text(size=fsize)) +
	annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) +
	annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1) +
	guides(fill=FALSE) +
	scale_fill_manual(values=c("indianred3", "royalblue3", "orange2")) 
	return(p1)
}


setwd("path/to/working/directory/") # set working directory

## Set helpful variables
metadata_file <- "[metadata_file].txt"
corrected_data_file <- "corrected_counts/gene/blood/[gene_counts].txt" # corrected data file
medz_ind_filt <- 100 # thresholds for defining global outlier
global_zscore <- 4 # z score for defining outlier for global outlier count

## Read in data 
gene_zscore <- fread(corrected_data_file, sep=",", header=TRUE, stringsAsFactor=FALSE)

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

## Set **percentile** threshold iterator
iterator <- c(0.5, 0.25, 0.1, (10/nrow(ids_filter)), (5/nrow(ids_filter)), (2/nrow(ids_filter)), (1/nrow(ids_filter))) # percentile thresholds
iterator <- unique(c(iterator, 1-iterator))
iterator <- iterator[order(iterator)]
iterator_length <- length(iterator)

outlier_indices <- list() # Get expression outliers
for (i in 1:length(iterator)) outlier_indices[[i]] <- pick_outliers(gene_zscore_scale_filter, iterator[i])

## Start up a parallel cluster
parallelCluster <- parallel::makeCluster(parallel::detectCores()) # adjust number of cores based on availability 

## Define worker function
wrapWorker <- function(gene_zscore_scale_filter, ids_filter, constraint_metrics, iterator, outlier_indices) {
	## Make data available to function
	force(gene_zscore_scale_filter)
	force(ids_filter)
	force(constraint_metrics)
	force(iterator)
	force(outlier_indices)

	# Count outliers for each case and control
	count_outlier_ids <- function(outlier_row, ids_filter, cases, controls) {
		outliers_in_gene <- ids_filter$sample_id[as.numeric(unlist(outlier_row))] # sample IDs of outliers
		sum_case <- length(which(outliers_in_gene %in% cases))
		sum_control <- length(which(outliers_in_gene %in% controls))
		return(c(sum_case, sum_control, sum_case / (sum_case + sum_control)))
	}

	## Find case and control outliers and test in logistic regression model
	outlier_enrichment <- function(iterator, ids_filter, gene_zscore_scale_filter, 
										constraint_metrics, random_labels_list, outlier_indices) {
		## Randomly reassign case/control
		ids_filter$affected_status <- random_labels_list

		case <- ids_filter$sample_id[ids_filter$affected_status == "Case"]
		control <- ids_filter$sample_id[ids_filter$affected_status == "Control"]

		## Initialize list
		pred <- c("lof_z", "mis_z", "syn_z")
		coeff_list <- list()
		for (p in pred) coeff_list[[p]] <- rep(NA, length(iterator))

		for (k in 1:length(iterator)) {
			## Get case/control outlier counts for each gene
			outlier_counts <- lapply(outlier_indices[[k]], function(outlier_current) count_outlier_ids(outlier_current, ids_filter, case, control))
			outlier_counts_matrix <- matrix(unlist(outlier_counts), nrow=length(outlier_counts), ncol=3, byrow=T)
			outlier_counts_df <- data.frame(outlier_counts_matrix)
			outlier_counts_df$GeneID <- names(outlier_indices[[k]])
			colnames(outlier_counts_df) <- c("CaseOutliers", "ControlOutliers", "ProportionCase", "GeneID")

			## Combine with constraint metrics
			constraint_metrics_merge <- merge(outlier_counts_df, constraint_metrics, by.x="GeneID", by.y="ensgene", all=FALSE)
			constraint_metrics_merge <- constraint_metrics_merge[, -1]

			## Remove zero counts (no outliers)
			constraint_metrics_merge <- constraint_metrics_merge[apply(constraint_metrics_merge, 1, function(row) sum(row[1:2])>0), ]

			## Logistic regression
			for (l in pred) {
				frm <- as.formula(paste("cbind(CaseOutliers, ControlOutliers) ~", l))
				model <- glm(formula=frm, data=constraint_metrics_merge, family=binomial)
				mod_summary <- summary(model)
				coeff_list[[l]][k] <- mod_summary[[12]][2, 1]
			}
		}
		return(coeff_list)
	}

	## Return results
	worker <- function(random_labels_list) {
    	outlier_enrichment(iterator, ids_filter, gene_zscore_scale_filter, 
    						constraint_metrics, random_labels_list, outlier_indices)
 	}
  	return(worker)
}

## Start parallel processing
npermute <- 10000 # number of permutations
random_labels <- replicate(npermute, sample(ids_filter$affected_status, nrow(ids_filter), replace=F))
random_labels_list <- split(random_labels, rep(1:ncol(random_labels), each=nrow(random_labels)))
outlier_permute <- parallel::parLapply(parallelCluster, random_labels_list, wrapWorker(gene_zscore_scale_filter, ids_filter, constraint_metrics, iterator, outlier_indices))

## Get real coefficients
real_coeff <- parallel::parLapply(parallelCluster, list(ids_filter$affected_status), wrapWorker(gene_zscore_scale_filter, ids_filter, constraint_metrics, iterator, outlier_indices))

## Write real coefficients in helpful format
pred <- c("lof_z", "mis_z", "syn_z")
for (k in pred) {
	unlist_coeff <- lapply(real_coeff, function(m) m[[k]])
	assign(paste0("coeff_", k), as.numeric(matrix(melt(unlist_coeff)[, 1], ncol=1, nrow=iterator_length, byrow=T)))
}

## Unlist permutations and compare with real coefficients
plot_collect <- data.frame(percentile=integer(),
							mean_coefficient=double(),
							error_low=double(),
							error_high=double(),
							emp_p=double(),
							real_coefficient=double(),
							predictor=character(),
							stringsAsFactors=FALSE)
	
for (i in pred) {
	tmp <- data.frame(array(NA, c(length(iterator), 9)))
	colnames(tmp) <- c("threshold", "mean_coefficient", "error_low", "error_high",
						 "emp_p", "real_coefficient", "predictor", "ngenes", "nsamples")

	unlist_coeff <- lapply(outlier_permute, function(m) m[[i]])
	matrix_coeff <- matrix(melt(unlist_coeff)[, 1], ncol=iterator_length, nrow=npermute, byrow=T)

	for (j in 1:iterator_length) {
		tmp[j, 1] <- j
		mean_coeff <- mean(matrix_coeff[, j], na.rm=T)
		sem <- sd(matrix_coeff[, j], na.rm=T) # empirical standard error
		tmp[j, 2] <- mean_coeff
		tmp[j, 3] <- mean_coeff-(1.96*sem)
		tmp[j, 4] <- mean_coeff+(1.96*sem)
		tmp[j, 5] <- (sum(abs(matrix_coeff[, j])>=abs(get(paste0("coeff_", i))[j]))+1) / (npermute+1) # compute empirical p-value
	}

	tmp$real_coefficient <- get(paste0("coeff_", i))
	tmp$predictor <- i
	tmp$nsamples <- nrow(ids_filter) # number of samples

	plot_collect <- rbind.data.frame(plot_collect, tmp)
}

## Plot
plot_collect$predictor <- revalue(plot_collect$predictor, c("lof_z"="Loss of Function", "mis_z"="Missense", "syn_z"="Synonymous"))
n <- plot_collect$nsamples[1]
emp <- emp_ast(plot_collect$emp_p)
label_data <- get_plot_labels(n)
ggsave("logistic_parallel.pdf", ggplot_build(plot_collect, round(label_data[2, ]), round(label_data[1, ], 3), emp), width=8.5, height=3.5)

## Shutdown cluster
if(!is.null(parallelCluster)) {
  parallel::stopCluster(parallelCluster)
  parallelCluster <- c()
}


