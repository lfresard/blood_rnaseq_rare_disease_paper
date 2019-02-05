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

# Master directory
dir = Sys.getenv('RARE_DIS_DIR')
## FUNCTIONS


load("figures_revision/Figure3.out.RData")

# Plot for solved case: splicing example


kctd7.df=res_splicing_number.df %>% filter(SAMPLE=="RD047")

kctd7.df$category=c("Rare variant", rep("Rare variant and filters",7),"HPO", "Splicing outlier", rep("Splicing and filters",5))
kctd7.df$category=factor(kctd7.df$category, levels=c("Splicing and filters","Splicing outlier","HPO","Rare variant and filters","Rare variant"))


kctd7.genes=get_genes_filters_splicing("RD047")

kctd7.genes=kctd7.genes[-1]

kctd7.df$is_causal_in=as.factor(unlist(lapply(kctd7.genes, is_causal,"KCTD7")))

#levels(kctd7.df$filter)=c("RV","RV_JUNC","RV_CADD","RV_HPO","RV_JUNC_CADD","RV_JUNC_HPO","RV_CADD_HPO","RV_JUNC_CADD_HPO","HPO","SPLI_OUTLIER", "SPLI_OUTLIER_PLI","SPLI_OUTLIER_HPO","SPLI_OUTLIER_RV", "SPLI_OUTLIER_RV_CADD","SPLI_OUTLIER_RV_CADD_HPO")

kctd7.dotchart=ggdotchart(kctd7.df, x = "filter", y = "value",color = "category",shape='is_causal_in',
	palette=c("#D0C6B1","#FFD800", "#587058", "#587498", "#E86850"),
	sorting = "descending",add = "segments",rotate = TRUE,
	group = "category",dot.size =20,label = round(kctd7.df.m$value),
	font.label = list(color = "white", size = 18, vjust = 0.5), ggtheme = theme_pubr(),
	add.params = list(color = "category", size = 2), 
	#legend=c(0.8,0.5), ylab="Number of candidate genes", xlab="Filter")#+RD_theme
	legend="", ylab="Number of candidate genes", xlab="Filter")#+RD_theme
kctd7.dotchart
KCTD7_rank.df <- data.frame(filter=c("Splicing outliers - Intitial count", "pLI over 0.9","HPO term matches gene", "Rare variant within 20bp","RV + CADD score over 10","RV + CADD score over 10 + HPO term match"), rank=c(8,NA,5,2,1,1))



KCTD7_rank.df$filter=factor(KCTD7_rank.df$filter, levels=c( "RV + CADD score over 10 + HPO term match","RV + CADD score over 10","Rare variant within 20bp","HPO term matches gene","pLI over 0.9", "Splicing outliers - Intitial count"))

KCTD7_rank.plot=ggplot(KCTD7_rank.df, aes(x=filter, y=rank)) + scale_y_continuous(breaks = pretty(KCTD7_rank.df$rank, n = 11))+RD_theme+
  geom_point(col="#D0C6B1", size=20) +   # Draw points
  geom_segment(x=1:6, xend=1:6, y=1, yend=11, linetype="dashed", size=0.1) + labs(y="zscore rank", x="") +theme(axis.text.y=element_blank())+ # Draw dashed lines
  #labs(title="Dot Plot", subtitle="Make Vs Avg. Mileage", caption="source: mpg") +  
  coord_flip()

pdf("/srv/scratch/restricted/rare_diseases/analysis/manuscript/figures_revision/kctd7.rank.pdf", height=5, width=15)
KCTD7_rank.plot
dev.off()
## combined_plot
kctd7_big_plot=ggdraw()+draw_plot(kctd7.dotchart, 0,1/3,2/3,2/3)+
	draw_plot(KCTD7_rank.plot, 2/3,1/3,1/3,1/3)+
	draw_plot_label(c('A', 'B', 'C'), c(0,0,1/2), c(1,1/3,1/3), size = 20)

pdf('/srv/scratch/restricted/rare_diseases/analysis/manuscript/figures_revision/kctd7_big_plot.pdf', w=10, h=12)
kctd7_big_plot
dev.off()

# Save data
save.image(file = paste(dir,"/analysis/manuscript/figures_revision/Figure6.in.RData",sep=""))
