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

## Make plots

# S4A
figS4a_dat <- read.table("[figS4a_dat].txt", header=T)
figS4a_plot <- ggplot(figS4a_dat, aes(x=sv, y=residuals)) + 
	geom_hline(yintercept=0) +
	geom_hline(yintercept=2, colour="gray", linetype="dashed") +
	geom_hline(yintercept=-2, colour="gray", linetype="dashed") +
	geom_point(aes(fill=batch), size=1.5, shape=21, colour="black") +
	facet_wrap( ~ cohort, scales="fixed") +
	labs(x="SV2", y="Residuals") +
	RD_theme +
	theme(panel.spacing.x=unit(1, "lines")) +
	annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) +
	annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1) +
	scale_fill_brewer(palette="Set3")


# S4B
# *TO FOLLOW*

# S4C
figS4c_dat <- read.table("[figS4c_dat].txt", header=T)
figS4c_plot <- ggplot(figS4c_dat, aes(pvalue)) +
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

## Write combined plots
sup_plot=ggdraw() + draw_plot(figS4a_plot, 0, 0.5, 1, 0.5) +
	draw_plot(figS4b_plot, 0, 0, 0.5, 0.5) +
	draw_plot(figS4c_plot, 0.5, 0, 0.5, 0.5) +
	draw_plot_label(c('A', 'B', 'C'), c(0, 0, 0.5), c(1, 0.5, 0.5), size=15)

pdf("sup_splines.pdf", w=7, h=6)
sup_plot
dev.off()


