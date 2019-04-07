#!/bin/R

## Code to generate Figure 2

library(cowplot)
library(ggplot2)
library(ggbeeswarm)

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

# function returning quartiles
median.quartile <- function(x){
	out <- quantile(x, probs = c(0.25,0.5,0.75))
	names(out) <- c("ymin","y","ymax")
	return(out) 
}

# Load data
fig2a=read.table("fig2a_data.txt", header=T, sep="\t")
fig2b=read.table("fig2b_data.txt", header=T, sep="\t")
fig2c=read.table("fig2c_data.txt", header=T, sep="\t")
fig2d=read.table("fig2d_data.txt", header=T, sep="\t")


## Make plots

# 2a

fig2a_plot <- ggplot(fig2a, aes(x=quantile, y=coefficient)) + geom_hline(yintercept = 0) + 
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


# 2b

fig2b_plot <- ggplot(fig2b) + 
	geom_boxplot(aes(x=Cutoff, y=OutlierCount+1, group=interaction(Cutoff, Direction), colour=Direction),
	 width=0.8, position=position_dodge(width=0.85)) +
	labs(x="Z-score", y="Number of Outlier Genes") +
	RD_theme +
	scale_y_log10(breaks=c(0,1,10,100,500,2000)) +
	scale_colour_manual(values=c("gray40", "gray1"))

# 2c
fig2c_plot=ggplot(fig2c,aes(x=Var1, y=value, fill=Var1))+
	geom_boxplot(color="black", notch=F, show.legend = FALSE)+
	geom_point(size = 0, stroke = 0)+
	scale_fill_manual(values=rep("lightgrey", 7), breaks=c("EXP_OUTLIER_PLI","EXP_OUTLIER_ASE","EXP_OUTLIER_RIVER", "EXP_OUTLIER_HPO","EXP_OUTLIER_RV", "EXP_OUTLIER_RV_CADD", "EXP_OUTLIER_RV_CADD_HPO"),
labels=c(expression("pLI " >="0.9"), "ASE", expression("RIVER " >= "0.85"), "HPO match", expression("Rare variant within 10kb","5 + CADD score  " >= "10"),"4 + 6"), name="Filter")+
	scale_x_discrete(breaks=c("EXP_OUTLIER_PLI","EXP_OUTLIER_ASE","EXP_OUTLIER_RIVER", "EXP_OUTLIER_HPO","EXP_OUTLIER_RV", "EXP_OUTLIER_RV_CADD", "EXP_OUTLIER_RV_CADD_HPO"),
labels=c("1", "2", "3","4","5", "6", "7"))+
	labs(x="Filter", y="Proportion of outlier genes") + RD_theme +guides(fill = guide_legend(override.aes = list(size = 4,shape =c(49,50,51,52,53,54,55))))+theme(legend.position = c(0.8,0.7))



# 2d

color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
mycol=sample(color, 75)

fig2d_plot= ggplot(fig2d, aes(filter,value))+ 
 	geom_quasirandom(aes(y =value, x = filter, colour = sample))+
 	stat_summary(fun.data = median.quartile, geom = "pointrange",position=position_nudge(x=0.5,y=0))+
	coord_flip()+
	#geom_path(group=sample, colour="grey")+
 	theme(legend.position="")+
 	scale_colour_manual(values=mycol) +
 	labs(y="log10(Number of genes +1)",x="Filter")+
 	scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) 

## Write combined plots
fig2_plot=ggdraw() + draw_plot(fig2a_plot, 0, 0.64, 1, 0.36) +
	draw_plot(fig2b_plot, 0, 0.39, 1, 0.25) +
	draw_plot(fig2c_plot, 0, 0, 0.5, 0.37) +
	draw_plot(fig2d_plot, 0.5, 0, 0.5, 0.37) +
	draw_plot_label(c('a', 'b', 'c', 'd'), c(0, 0, 0, 0.5), c(1, 0.63, 0.38, 0.38), size=15)

pdf("fig2.pdf", w=8.5, h=10)
fig2_plot
dev.off()


