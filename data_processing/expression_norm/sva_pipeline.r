#!/bin/R

# sva.r
# Description: Performs Surrogate Variable Analysis (SVA) on matrix of RNA-Seq count data output by RSEM
# Outputs: PCA plot of corrected and uncorrected data, corrected matrix, matrix of significant surrogate variables (SVs), QC plots

## Load required libraries
library(data.table)
library(sva)
library(ggplot2)
library(ggfortify)
library(tximport)
library(dplyr)
library(readr)
library(stringr)
library(splines)
library(RColorBrewer)
library(plyr)

## Function to return all SVs associated with a given predictor, at specified p-value cutoff
predictor_sig <- function(sva_object, predictor, cutoff) {
	return(which(apply(sva_object, 2, function(x) 
		any(summary(lm(x ~ factor(predictor)))[[4]][-1,4] < cutoff)))) 
}

## Set working directory
setwd("/path/to/working/directory/")

## Load sample info
metadata_file <- "[metadata_file].txt"
metadata <- read_tsv(metadata_file) %>% filter(in_freeze=="yes") %>% filter(source_of_RNA=="Blood")
list_DGN_counts <- list.files(path="/path/to/rsem/gene/counts/", pattern=glob2rx("LD*.genes.results*"))
DGN_names <- str_extract(list_DGN_counts, "LD[0-9]+")
sample_id <- c(metadata$sample_id, DGN_names)

# Get list of files from directory containing RSEM output
# assumes one *.genes.results file per sample
dir <- "/path/to/rsem/gene/counts/"
file_list <- list.files(dir)

## Read data
# Read first file
index <- file_list[grep(paste0(sample_id[1], ".genes.results"), file_list)]
dat <- fread(paste0(dir, index), header=T, sep="\t", stringsAsFactor=F)
dat <- dat[ ,c("gene_id","TPM")] # keep gene name, TPM
colnames(dat) <- c("gene_id", sample_id[1])

# Read remaining files
for (i in 2:length(sample_id)) {
	cat(paste0(i, " "))
	index_temp <- file_list[grep(paste0(sample_id[i], ".genes.results"), file_list)]
	if (length(index_temp) != 0) {
		temp <- fread(paste0(dir, index_temp), header=T, sep="\t", stringsAsFactor=F)
		stopifnot(temp$Geneid==dat$gene_id) # check rownames match, throw error if not
		dat <- cbind(dat, temp[, "TPM"]) # add column to data matrix
		colnames(dat)[ncol(dat)] <- sample_id[i] # change name of appended column
	}
}

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

## SVA
mod <- model.matrix(~1, data=as.data.frame(t(dat_filter_log_scale)))
sva_fit <- sva(dat_filter_log_scale, mod, method="two-step")

## Splines
batch <- c(metadata$batch, rep(99, length(DGN_names)))
study <- c(metadata$institution, rep("DGN", length(DGN_names)))
sig_sv <- unique(c(predictor_sig(sva_fit$sv, batch, 1e-30), predictor_sig(sva_fit$sv, study, 1e-30)))

## Regress out SVs
modsv <- cbind(mod, sva_fit$sv)
modsv_no_spline <- modsv
for (i in sig_sv) modsv <- cbind(modsv, bs(sva_fit$sv[, i], df=60, degree=1))
fitsv <- lm.fit(modsv, t(dat_filter_log_scale))

## Write corrected data
write_dat <- cbind(rownames(fitsv$residuals), fitsv$residuals)
colnames(write_dat)[1] <- "gene_name"
write.table(write_dat, "gene_counts_corrected_freeze_spline.txt", sep=",", row.names=F, col.names=T) # write corrected matrix
write.table(as.data.frame(sva_fit[1]), "gene_sig_svs_freeze_spline.txt", sep=",", row.names=F, col.names=T) # write surrogate variables

## Likehood ratio test
lrt_fit <- apply(t(dat_filter_log_scale), 2, function(x) {
	collect <- matrix(NA, nrow=1, ncol=5)
	sv <- lm(x~modsv_no_spline)
	sv_spline <- lm(x~modsv)
	collect[1] <- anova(sv, sv_spline, test="LRT")[[5]][2]
	sv_summary <- summary(sv)
	sv_spline_summary <- summary(sv_spline)
	collect[2:3] <- c(sv_summary$r.squared, sv_summary$adj.r.squared)
	collect[4:5] <- c(sv_spline_summary$r.squared, sv_spline_summary$adj.r.squared)
	return(collect)})

lrt_fit <- data.frame(t(lrt_fit))
colnames(lrt_fit) <- c("pvalue", "r2_no_spline", "adj_r2_no_spline", "r2_spline", "adj_r2_spline")
lrt_fit$r2_diff <- lrt_fit$r2_spline - lrt_fit$r2_no_spline
lrt_fit$adj_r2_diff <- lrt_fit$adj_r2_spline - lrt_fit$adj_r2_no_spline

adj_pvals <- p.adjust(lrt_fit[,1], method="BH")
adj_pvals_thresh <- sum(adj_pvals<0.05)

## Plot LRT pvalue result
lrt_plot <- data.frame(factor(1:nrow(lrt_fit)), lrt_fit[, 1])
colnames(lrt_plot) <- c("gene", "pvalue")

p_lrt <- ggplot(lrt_plot, aes(pvalue)) +
	geom_histogram(binwidth=0.02, fill="gray40") +
	labs(x="P Value", y="Count") +
	theme_bw() +
	theme(strip.background=element_blank(),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),
		axis.text=element_text(size=9),
		panel.border=element_blank()) +
	annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) +
	annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1)
ggsave("lrt_plot.pdf", p_lrt, width=3.5, height=2)

## Plot R squared
r2_plot <- cbind(1:nrow(lrt_fit), lrt_fit[, c("r2_diff", "adj_r2_diff")])
colnames(r2_plot)[1] <- "gene_id"
r2_plot <- melt(r2_plot, id="gene_id")

plot_current <- subset(r2_plot, variable=="r2_diff")[order(subset(r2_plot, variable=="r2_diff")$value, decreasing=T), ]
plot_current$gene_id <- 1:nrow(plot_current)

p_r2 <- ggplot(data=plot_current, aes(x=gene_id, y=value)) + #
	geom_bar(stat="identity", fill="blue") +
	geom_hline(yintercept=0) +
	labs(x="Genes", y="Change in R^2") +
	scale_x_continuous(breaks=c(1, seq(5000, nrow(lrt_fit), 5000))) +
	theme_bw() +
	theme(strip.background=element_blank(),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),
		axis.text=element_text(size=9),
		panel.border=element_blank()) +
	annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) +
	annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1) 
ggsave("r2_plot.pdf", p_r2, width=4, height=3)

## Plot PCA
pca_object <- prcomp(fitsv$residuals, center=TRUE, scale=TRUE)
combined_nonnum <- cbind(fitsv$residuals, as.data.frame(factor(batch))) # non-numeric holder
colnames(combined_nonnum)[ncol(combined_nonnum)] <- "batch"

p1 <- autoplot(pca_object, scale=0, data=combined_nonnum, colour="batch") +
	geom_point(aes(colour=batch), alpha=c(rep(1, length(grep("RD", rownames(combined_nonnum)))), rep(0, length(grep("LD", rownames(combined_nonnum)))))) +
	theme_bw() +
	theme(strip.background=element_blank(),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),
		axis.text=element_text(size=10),
		panel.border=element_blank()) +
	annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) +
	annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1) +
	guides(colour=FALSE) +
	scale_colour_brewer(palette="Set3") + xlim(-310, 200) + ylim(-230, 60) # optionally manually set axes length
ggsave("gene_sva_freeze.pdf", p1, width=3, height=3)

## Plot PCA (no correction)
pca_object_none <- prcomp(t(dat_filter_log_scale), center=TRUE, scale=TRUE)

p2 <- autoplot(pca_object_none, scale=0, data=combined_nonnum, colour="batch") +
	theme_bw() +
	theme(strip.background=element_blank(),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),
		axis.text=element_text(size=10),
		panel.border=element_blank()) +
	annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) +
	annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1) +
	guides(colour=FALSE) +
	scale_colour_brewer(palette="Set3") + xlim(-310, 200) + ylim(-230, 60) # optionally manually set axes length
ggsave("gene_no_sva_freeze.pdf", p2, width=3, height=3)

## Load DGN metadata
DGN_data <- load("/path/to/DGN/metadata/")
meta_DGN <- combinedEnv
meta_DGN <- subset(meta_DGN, rownames(meta_DGN) %in% DGN_names)
meta_DGN_df <- data.frame(meta_DGN)

covariates=data.frame(institution=as.integer(as.factor(study)),
	batch=as.integer(as.factor(batch)),
	age=c(metadata$age, meta_DGN_df$Agegroup), 
	sex=as.integer(as.factor(c(ifelse(as.integer(as.factor(metadata$sex))==1, 2, 1), meta_DGN_df$Sex))), 
	affected_status=as.numeric(as.factor(c(metadata$affected_status, rep("Control", length(DGN_names))))),
	read_length=c(as.integer(as.factor(metadata$read_length)), rep(4, length(DGN_names))),
	sequencer=c(as.integer(as.factor(c(rep("nextseq", nrow(metadata)), rep("highseq", length(DGN_names)))))),
	read_type=c(as.integer(as.factor(c(rep("paired", nrow(metadata)), rep("single", length(DGN_names)))))))

# Split batch
batch=as.factor(as.integer(as.factor(batch)))
batch_mat=model.matrix(~batch)
batch_mat=as.data.frame(batch_mat[, 2:ncol(batch_mat)])
batch_mat$batch1=c(ifelse(metadata$batch==1,1,0), rep(0, length(DGN_names)))

# Split study
study=as.factor(as.integer(as.factor(study)))
study_mat=model.matrix(~study)
study_mat=as.data.frame(study_mat[, 2:ncol(study_mat)])
study_mat$study1=ifelse(study==1,1,0)

# Split read length
read_length=as.factor(as.factor(c(as.integer(as.factor(metadata$read_length)), rep(4, length(DGN_names)))))
read_length_mat=model.matrix(~read_length)
read_length_mat=as.data.frame(read_length_mat[, 2:ncol(read_length_mat)])
read_length_mat$read_length1=ifelse(read_length==1,1,0)

covariates=cbind(covariates, batch_mat, study_mat, read_length_mat)
cor.pc.cov=cor(sva_fit$sv,covariates, use="complete.obs")

cor.pc.cov.m=melt(cor.pc.cov)
cor.pc.cov.m$Var2=factor(cor.pc.cov.m$Var2, levels=c("institution", "study1", "study2", "study3", "study4", "batch", "batch1", "batch2", "batch3", "batch4", "batch5", "batch6", "batch7", "age", "sex", "affected_status", "read_length", "read_length4", "read_length1", "read_length2", "read_type", "sequencer"))
cor.pc.cov.m <- subset(cor.pc.cov.m, Var2!="batch")
cor.pc.cov.m <- subset(cor.pc.cov.m, Var2!="institution")
cor.pc.cov.m <- subset(cor.pc.cov.m, Var2!="read_length")
cor.pc.cov.m$Var2 <- revalue(cor.pc.cov.m$Var2, c("study1"="CGS", "study2"="CHEO", "study3"="DGN", "study4"="UDN"))
cor.pc.cov.m$Var2 <- revalue(cor.pc.cov.m$Var2, c("read_length1"="read_75bp", "read_length2"="read_150bp", "read_length4"="read_50bp"))

# Plot heatmap
cor.heatmap=ggplot(data=cor.pc.cov.m, aes(x=Var2, y=as.factor(Var1), fill=value)) + 
 	geom_tile(color="white") +
 	scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0, limit=c(-1,1), space="Lab", name="Pearson\nCorrelation") +
 	theme_minimal() + 
 	coord_fixed() + labs(x="Covariates", y="Surrogate Variables") + 
 	theme(legend.position="right",
 		axis.text.x=element_text(angle=45, hjust=1),
 		axis.text=element_text(size=9)) +
 	annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) +
	annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1)
ggsave("sup_exp_corPCcov.pdf", cor.heatmap, width=5, height=5)

## Combined_plot
sup_pca_plot=ggdraw() + draw_plot(p2, 0, 0.6, 0.5, 0.4) +
	draw_plot(p1, 0.5, 0.6, 0.5, 0.4) +
	draw_plot(cor.heatmap, 0, 0, 1, 0.6) +
	draw_plot_label(c('A', 'B', 'C'), c(0, 0.5, 0), c(1, 1, 0.6), size=15)

pdf("data_correction_sup.pdf", w=7, h=8)
	sup_pca_plot
dev.off()

