#!/bin/R

## Code to generate Figure S4

library(cowplot)
library(ggplot2)

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

## Function to return numeric p-value in asterick notation
emp_ast <- function(pvalues) {
	sapply(pvalues, function(i) {
		if (i > 0.05) return("")
		else if (i <= 0.5 & i > 0.01) return("*")
		else if (i <= 0.01 & i > 0.001) return("**")
		else return("***")
		})
}

## Make plots

# S5A
# **Overview image imported from external software. Space is reserved when writing combined plot**

# S5B
figS5b_dat <- read.table("[figS5b_dat].txt", header=T)
figS5b_plot <- ggplot(figS5b_dat, aes(x=factor(percentile), y=value)) +
	geom_hline(yintercept=0) +
	geom_hline(yintercept=2, colour="gray", linetype="dashed") +
	geom_hline(yintercept=-2, colour="gray", linetype="dashed") +
	scale_y_continuous(breaks=seq(-10, 10, 2), labels=seq(-10, 10, 2)) +
	labs(x="Percentile", y="Z-score") +
	geom_boxplot(colour="gray50") +
	RD_theme +
	theme(axis.text.x=element_text(angle=45, hjust=1)) +
	annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) +
	annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1)

# S5C
figS5c_dat <- read.table("[figS5c_dat].txt", header=T)
figS5c_plot <- ggplot(figS5c_dat, aes(x=quantile, y=coefficient)) + geom_hline(yintercept=0) + 
	geom_pointrange(aes(ymin=error_low, ymax=error_high, colour=predictor)) +
	facet_grid(. ~ predictor, scales="fixed") +
	labs(x="Percentile", y="Log Odds") +
	geom_text(aes(x=quantile, y=error_high+0.004, label=emp_ast(pvalue)), size=3) +
	RD_theme +
	theme(axis.text.x=element_text(angle=45, hjust=1),
		strip.text=element_text(size=fsize)) +
	annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) +
	annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1) +
	guides(colour=FALSE) +
	scale_colour_manual(values=c("indianred3", "royalblue3", "orange2"))

# S5D
figS5d_dat <- read.table("[figS5d_dat].txt", header=T)
figS5d_plot <- ggplot(figS5d_dat) + geom_hline(yintercept=0) + 
	geom_errorbar(aes(x=threshold, ymin=error_low, ymax=error_high), colour="gray") + 
	geom_point(aes(x=threshold, y=real_coefficient, fill=predictor), shape=21, colour="Black") +
	facet_grid(. ~ predictor, scales="free_y") +
	labs(x="Percentile", y="Log Odds") +
	geom_text(aes(x=threshold, y=real_coefficient+0.01, label=emp_ast(pvalue), size=3) +
	RD_theme +
	theme(axis.text.x=element_text(angle=45, hjust=1),
		strip.text=element_text(size=fsize)) +
	annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) +
	annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1) +
	guides(fill=FALSE) +
	scale_fill_manual(values=c("indianred3", "royalblue3", "orange2")) 

# S5E
figS5e_dat <- read.table("[figS5e_dat].txt", header=T)
figS5e_plot <- ggplot(figS5e_dat, aes(x=quantile, y=coefficient)) + geom_hline(yintercept=0) + 
	geom_pointrange(aes(ymin=error_low, ymax=error_high, colour=predictor)) +
	facet_grid(. ~ predictor, scales="fixed") +
	labs(x="Z-score", y="Log Odds") +
	geom_text(aes(x=quantile, y=error_high+0.004, label=emp_ast(pvalue)), size=3) +
	RD_theme +
	theme(axis.text.x=element_text(angle=45, hjust=1),
		strip.text=element_text(size=fsize)) +
	annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) +
	annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1) +
	guides(colour=FALSE) +
	scale_colour_manual(values=c("indianred3", "royalblue3", "orange2"))


## Write combined plots
sup_plot=ggdraw() + 
	draw_plot(figS5b_plot, 0.5, 0.74, 0.5, 0.26) +
	draw_plot(figS5c_plot, 0, 0.47, 1, 0.27) +
	draw_plot(figS5d_plot, 0, 0.2, 1, 0.27) +
	draw_plot(figS5e_plot, 0, 0, 1, 0.2) +
	draw_plot_label(c('A', 'B', 'C', 'D', 'E'), c(0, 0.49, 0, 0, 0), c(1, 1, 0.74, 0.47, 0.2), size=15)

pdf("sup_logistic_model.pdf", w=16, h=15)
sup_plot
dev.off()
