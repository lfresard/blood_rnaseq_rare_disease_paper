#!/bin/R

## Runs one iteration of logistic regression model with varying number of controls as shown in Sup Fig 6

## Load required libraries
library(data.table)
library(sva)
library(ggplot2)
library(ggfortify)
library(splines)
library(annotables)

## <------------------- FUNCTIONS

## Function to get surrogate variables significantly associated with a given predictor
predictor_sig <- function(sva_object, predictor, cutoff) {
        return(which(apply(sva_object, 2, function(x)
                any(summary(lm(x ~ factor(predictor)))[[4]][-1,4] < cutoff))))
}

## Function to get all outliers at a specified percentile threshold
pick_outliers <- function(medz, quantile_thresh, method) {
        if (method == "<") return(apply(t(medz), 1, function(for_row) {
                        get_quantile <- as.numeric(quantile(for_row, quantile_thresh))
                        list(which(for_row <= get_quantile))
                 }))
        else if (method == ">") return(apply(t(medz), 1, function(for_row) {
                        get_quantile <- as.numeric(quantile(for_row, quantile_thresh))
                        list(which(for_row >= get_quantile))
                 }))
        else print("Error: choose '<' or '>'")
}

## Function to count number of outliers for cases and controls
count_outlier_ids <- function(outlier_row, ids_filter, cases, controls) {
        outliers_in_gene <- ids_filter$sample_id[as.numeric(unlist(outlier_row))] # sample IDs of outliers
        sum_case <- length(which(outliers_in_gene %in% cases))
        sum_control <- length(which(outliers_in_gene %in% controls))
        return(c(sum_case, sum_control, sum_case / (sum_case + sum_control)))
}

## <------------------- MAIN

setwd("[path/to/working/directory/]") # set working directory

## Process input parameter
args <- commandArgs(trailingOnly=TRUE)
random_indices_vector <- as.numeric(args[1]) # input is random index column to choose for current run

## Load sample metadata
ids <- read.table("[metadata_file].txt", sep=",", header=T, stringsAsFactor=F)

# Load raw gene count data
dat <- fread("[raw_gene_counts].txt", sep=",", header=T)

## Load random indices for data subset
dgn_random <- read.table("[random_control_samples_indices].txt", sep="\t", header=FALSE)

global_out <- [global_out].txt # file containing global outlier sample IDs

## Copy data file to restore from after random subset
dat_copy <- dat
sample_id_copy <- ids

n_dgn_control <- c(100, 200, 300, 400, 500, 750, 900)
df_spline <- c(11, 17, 24, 30, 36, 51, 60) # degrees of freedom for regression splines for each control N

coeff_list <- list()

for (p in 1:length(n_dgn_control)) {
        ## ------------------- START SURROGATE VARIABLE ANALYSIS
        dat <- dat_copy
        sample_id <- sample_id_copy

        ## Subset DGN control samples
        rd_index <- grep("RD", colnames(dat))
        dgn_sample <- grep("LD", colnames(dat))[dgn_random[1:n_dgn_control[p], random_indices_vector]]
        dat <- dat[ , c(1, rd_index, dgn_sample), with=F]
        sample_id <- sample_id$sample_id[which(sample_id$sample_id %in% colnames(dat))]

        ## Filter lowly-expressed genes in each sample origin (RD or LD)
        dat_filter_genes <- sapply(unique(substr(sample_id, 1, 2)), function(i) {
                        tmp <- dat[, c(TRUE, colnames(dat[, -1]) %in% sample_id[grep(i, sample_id)]), with=FALSE]
                        as.character(unlist(tmp[apply(tmp[, -1], 1, function(x) sum(x>0.5) > (length(colnames(tmp))*0.5)), 1]))
                })

        if (length(unique(substr(sample_id, 1, 2))) > 1) {
                dat_filter_genes_intersect <- Reduce(intersect, dat_filter_genes) # find common set of expressed genes
                } else dat_filter_genes_intersect <- dat_filter_genes

        dat_filter <- dat[which(unlist(dat[,1]) %in% dat_filter_genes_intersect), ]

        ## Moderated log transform
        dat_filter_log <- dat_filter[, apply(dat_filter[, -1], 2, function(x) log10(x+1))]
        dat_filter_log <- cbind(dat_filter[, 1], dat_filter_log)

        ## Unit variance and center for cols with var != 0
        temp_dat <- t(dat_filter_log[ ,-1])
        colnames(temp_dat) <- as.character(unlist(dat_filter_log[,1]))

        for (i in 1:ncol(temp_dat)) {
                if (var(temp_dat[ ,i]) != 0) temp_dat[, i] <- scale(temp_dat[, i], center=TRUE, scale=TRUE)
        }

        temp_dat_var <- temp_dat[, apply(temp_dat, 2, function(x) !var(x) == 0)] # remove genes with zero variance

        dat_filter_log_scale <- t(temp_dat_var) # rename and transpose

        ## Perform SVA
        mod <- model.matrix(~1, data=as.data.frame(t(dat_filter_log_scale)))
        sva_fit <- sva(dat_filter_log_scale, mod, method="two-step")

        ## Regression splines
        batch <- c(ids$batch[grep("RD", ids$sample_id)], rep(10, length(dgn_sample)))
        study <- c(ids$institution[grep("RD", ids$sample_id)], rep("DGN", length(dgn_sample)))
        sig_sv <- unique(c(predictor_sig(sva_fit$sv, batch, 1e-30), predictor_sig(sva_fit$sv, study, 1e-30)))

        ## Regress out SVs
        modsv <- cbind(mod, sva_fit$sv)
        for (i in sig_sv) modsv <- cbind(modsv, bs(sva_fit$sv[, i], df=df_spline[p], degree=1))
        fitsv <- lm.fit(modsv, t(dat_filter_log_scale))

        ## -------------------------- START LOGISTIC REGRESSION MODEL
        ## Scale and center matrix
        gene_zscore_scale <- scale(fitsv$residuals, center=TRUE, scale=TRUE)
        colnames(gene_zscore_scale) <- unlist(strsplit(as.character(colnames(gene_zscore_scale)), '[.]'))[c(TRUE,FALSE)]

        ## Save gene names to vector
        genes <- colnames(gene_zscore_scale)

        ## Get constraint metrics for current gene set
        constraint_metrics_all <- read.table("[ExAC_constraint_metrics].txt", sep="\t", header=TRUE, stringsAsFactor=FALSE)
        constraint_metrics <- merge(constraint_metrics_all[, c("gene", "syn_z", "mis_z", "lof_z", "pLI")], grch37[, c("ensgene", "symbol")], by.x="gene", by.y="symbol", all=FALSE)
        constraint_metrics <- constraint_metrics[-which(duplicated(constraint_metrics)), ] # drop duplicated rows
        constraint_metrics <- constraint_metrics[, c("ensgene", "pLI", "mis_z", "syn_z", "lof_z")] # re-arrange columns
        constraint_metrics <- constraint_metrics[which(constraint_metrics$ensgene %in% genes), ] # filter to include only genes in current corrected matrix

        for (s in c("no_switch", "switch")) {
                ## Remove global expression outliers from corrected data
                gene_zscore_scale_filter <- gene_zscore_scale[-which(rownames(gene_zscore_scale) %in% global_out), ] # remove global outliers from matrix
                ids_filter <- ids[which(ids$sample_id %in% rownames(gene_zscore_scale_filter)), ] # update IDs to remove IDs not in filtered matrix
                ids_filter <- ids_filter[match(rownames(gene_zscore_scale_filter), ids_filter$sample_id), ] # match IDs if necessary

                if (s=="switch") {
                        rd_index <- grep("Case", ids_filter$affected_status)
                        ids_filter$affected_status[grep("RD", ids_filter$sample_id)] <- "Case"
                        ids_filter$affected_status[rd_index] <- "Control"
                }

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
                plot_collect_df[, 1:7] <- apply(plot_collect_df[, c(1:7)], 2, function(for_col) as.numeric(for_col))
                plot_collect_df$n_dgn_control <- n_dgn_control[p]
                plot_collect_df$cohort <- s
                coeff_list[[length(coeff_list)+1]] <- plot_collect_df
        }

}

## Write output
coeff_unlist <- do.call("rbind", coeff_list)
coeff_unlist$run <- random_indices_vector
write.table(coeff_unlist, file=paste0("output/log_model_switch_run", random_indices_vector, ".txt"), sep="\t", col.names=T, row.names=F)

