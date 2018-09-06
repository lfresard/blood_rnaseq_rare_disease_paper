#!/bin/R

## Code to generate Figure 2

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

# 2A
2a_dat <- read.table("[figure_2a_data].txt", header=T)

fig2a <- ggplot(2a_dat, aes(x=quantile, y=coefficient)) + geom_hline(yintercept = 0) + 
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

# 2B
fig2b_dat <- read.table("[figure_2b_data].txt", header=F, sep="\t")
n_samples <- c("n=1", "n=2", "n=5", "n=10", "10%", "25%", "50%", "75%", "90%", "n=10", "n=5", "n=2", "n=1")
fig2b_dat_subset <- subset(log_test_data, predictor=="lof_z" & quantile<5)

p2 <- ggplot(fig2b_dat_subset, aes(x=factor(quantile), y=coefficient)) + 
	geom_hline(yintercept=0) + 
	geom_boxplot(aes(fill=cohort, alpha=factor(n_dgn_control)), colour="black", size=0.3, outlier.shape=1, outlier.size=0.5) +
	scale_x_discrete(labels=n_samples) + 
	labs(x="Percentile (Under Expression)", y="Log Odds", title="Loss of Function") +
	RD_theme +
	annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) +
	annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1) +
	guides(alpha=FALSE, colour=FALSE, fill=FALSE) +
	scale_fill_manual(values=c("indianred3"))

# 2C
fig2c_dat <- read.table("[figure_2c_data].txt", header=T)

fig2c <- ggplot(plot_collect) + 
	geom_boxplot(aes(x=Cutoff, y=OutlierCount+1, group=interaction(Cutoff, Direction), colour=Direction),
	 width=0.8, position=position_dodge(width=0.85)) +
	labs(x="Z-score", y="Number of Outlier Genes") +
	RD_theme +
	scale_y_log10(breaks=c(0,1,10,100,500,2000)) +
	scale_colour_manual(values=c("gray40", "gray1"))

# 2D
load("Figure2D.out.RData")
fig2d <- filter_withsglt.exp.case.10kb.plot + theme(legend.position="")

## Write combined plots
fig2_plot=ggdraw() + draw_plot(fig2a, 0, 0.64, 1, 0.36) +
	draw_plot(fig2b, 0, 0.39, 1, 0.25) +
	draw_plot(fig2c, 0, 0, 0.5, 0.37) +
	draw_plot(fig2d, 0.5, 0, 0.5, 0.37) +
	draw_plot_label(c('A', 'B', 'C', 'D'), c(0, 0, 0, 0.5), c(1, 0.63, 0.38, 0.38), size=15)

pdf("fig2.pdf", w=8.5, h=10)
fig2_plot
dev.off()


