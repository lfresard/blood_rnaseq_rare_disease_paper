#!/bin/R
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(reshape2)
library(data.table)
library(annotables)
library(ggpubr)
library(cowplot)
rm(list=ls())

is_causal=function(item,symbol){
	if (any(item==symbol)==TRUE){
		causal=1
	}else{
		causal=0
	}
	return(causal)
}


# Read in Figure 4 data
fig4a=read.table("Fig4a_data.txt",sep="\t", header=T)
# Plot for solved case: splicing example

fig4a_plot=ggdotchart(fig4a, x = "filter", y = "Number.of.candidates",color = "category",shape=as.factor('is_causal_in'),
	palette=c("#D0C6B1","#FFD800", "#587058", "#587498", "#E86850"),
	sorting = "descending",add = "segments",rotate = TRUE,
	group = "category",dot.size =20,label = round(fig4a$Number.of.candidates),
	font.label = list(color = "white", size = 18, vjust = 0.5), ggtheme = theme_pubr(),
	add.params = list(color = "category", size = 2), 
	#legend=c(0.8,0.5), ylab="Number of candidate genes", xlab="Filter")#+RD_theme
	legend="", ylab="Number of candidate genes", xlab="Filter")#+RD_theme


## Final plot was modified with Inkscape for aestetics purposes.