#!/bin/R

## Code to generate Figure S3

library(cowplot)
library(ggplot2)
library(RColorBrewer)

setwd("[path/to/working/directory/]") # set working directory

## Define ggplot theme
fsize <- 15
RD_theme = theme_classic() +
	theme(axis.text.x = element_text(size=fsize),
	axis.text.y = element_text(size=fsize), 
	axis.title = element_text(size=fsize), 
    legend.text = element_text(size=fsize-1), 
    legend.title = element_text(size=fsize), 
    axis.ticks = element_line(size=0.1))

## Make plots

# S3A
load("FigureS3A.out.RData")
figS3a_dat <- pca_corrected

fig_S3a <- autoplot(figS3a_dat, scale=0, colour="batch") +
	geom_point(aes(colour=batch), alpha=c(rep(1, length(grep("RD", rownames(figS3a_dat)))), rep(0, length(grep("LD", rownames(figS3a_dat)))))) +
	RD_theme +
	annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) +
	annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1) +
	guides(colour=FALSE) +
	scale_colour_brewer(palette="Set3")

# S3B
load("FigureS3B.out.RData")
figS3b <- pca_uncorrected

fig_S3b <- autoplot(figS3b_dat, scale=0, colour="batch") +
	geom_point(aes(colour=batch), alpha=c(rep(1, length(grep("RD", rownames(figS3b_dat)))), rep(0, length(grep("LD", rownames(figS3b_dat)))))) +
	RD_theme +
	annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) +
	annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1) +
	guides(colour=FALSE) +
	scale_colour_brewer(palette="Set3")

# S3C
figS3c_dat <- read.table("cor_sv_covar.txt", header=T)

fig_S3c <- ggplot(data=figS3c_dat, aes(x=Var2, y=as.factor(Var1), fill=value)) + 
 	geom_tile(color="white") +
 	scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0, limit=c(-1,1), space="Lab", name="Pearson\nCorrelation") +
 	theme_minimal() + 
 	coord_fixed() + labs(x="Covariates", y="Surrogate Variable") + 
 	theme(legend.position="right",
 		axis.text.x=element_text(angle=45, hjust=1),
 		axis.text=element_text(size=9)) +
 	annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) +
	annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1)

## Write combined plots
## Combined_plot
sup_pca_plot=ggdraw() + draw_plot(fig_S3a, 0, 0.6, 0.5, 0.4) +
	draw_plot(fig_S3b, 0.5, 0.6, 0.5, 0.4) +
	draw_plot(fig_S3c, 0, 0, 1, 0.6) +
	draw_plot_label(c('A', 'B', 'C'), c(0, 0.5, 0), c(1, 1, 0.6), size=15)

pdf("figS3.pdf", w=7, h=8)
	sup_pca_plot
dev.off()
