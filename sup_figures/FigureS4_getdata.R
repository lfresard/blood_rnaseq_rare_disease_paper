#!/bin/R

## Code to generate data for Figure S4

setwd("[path/to/working/directory/]") # set working directory

## Read and procress raw gene count data
gene_counts_raw <- read.table('[gene_count_raw_log_scale].txt', header=T) # raw counts of chosen gene to plot
ids <- read.table('[sample_metadata].txt', header=T)

gene <- data.frame(colnames(gene_counts_raw), "gene1", gene_counts_raw, 
	factor(ids$batch[match(colnames(dat_filter_log_scale), ids$sample_id)]), stringsAsFactors=F)
colnames(temp_gene) <- c("sample_id", "name", "exp", "batch")

## Get residual vector for no model, no splines, splines
temp1 <- temp_gene
temp1$residuals <- temp_gene$exp
temp1$cohort <- "NoModel"
temp1$sv <- sva_fit$sv[, 2]

temp2 <- temp_gene
temp2$residuals <- scale(lm(temp2$exp~cbind(mod, sva_fit$sv))$residuals, center=T, scale=T)
temp2$cohort <- "NoSplines"
temp2$sv <- sva_fit$sv[, 2]

temp3 <- temp_gene
mod_temp <- cbind(mod, sva_fit$sv, bs(sva_fit$sv[, 3], df=60, degree=1)) #cbind(mod, sva_fit$sv)
for (i in 1:length(sig_sv)) mod_temp <- cbind(mod_temp, bs(sva_fit$sv[, sig_sv[i]], df=60, degree=1)) 
temp3$residuals <- scale(lm(temp2$exp~modsv)$residuals, center=T, scale=T)
temp3$cohort <- "Splines"
temp3$sv <- sva_fit$sv[, 2]

temp_full <- rbind(temp1, temp2, temp3)

## Write data
write.table(temp_full, file="figS4a_dat.txt", row.names=F, col.names=T)