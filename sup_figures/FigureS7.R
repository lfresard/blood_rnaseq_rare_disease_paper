#!/bin/R

# plot the percent of sample(with gen data) for which there is at least one candidate accross filters



#Libraries
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(dplyr)
library(readr)
library(reshape2)

# Main
rm(list=ls())
load("/srv/scratch/restricted/rare_diseases/analysis/outlier_analysis/expression_level/expression_filters_genes.RData")
load("/srv/scratch/restricted/rare_diseases/analysis/outlier_analysis/splicing/for_freeze/splicing_filters_genes.RData")

fsize=20
RD_theme=theme_classic()+
	theme(axis.text.x= element_text(size=fsize),
	axis.text.y= element_text(size=fsize), 
	axis.title = element_text(size = fsize), 
    legend.text = element_text(size = fsize-1), 
    legend.title = element_text(size = fsize), 
    axis.ticks = element_line(size = 0.1),
	strip.background=element_blank(),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	axis.text=element_text(size=9),
	panel.border=element_blank()) 


# Get metadata
metadata = "/srv/scratch/restricted/rare_diseases/data/metadata/2018_06_12_Rare_Disease_Metadata.tsv"
metadata = read_tsv(metadata) %>% filter(in_freeze=="yes") %>% filter(source_of_RNA=="Blood")
sample_with_data=metadata %>% filter(variant_data=="exome" | variant_data=="genome")%>% select(sample_id) %>% pull



# Expression non zero individual in %
non_zero_samples=colSums(sample_exp_outlier_filter_withsgl_10kb_up.df[,3:(ncol(sample_exp_outlier_filter_withsgl_10kb_up.df)-3)] != 0)/length(sample_with_data) *100
non_zero_samples=melt(non_zero_samples)
non_zero_samples$filter=c(1:13)
 

exp_non_zero_samples_plot=ggplot(non_zero_samples, aes(x=as.factor(filter), y=value)) + 
  geom_point(size=3) + 
  geom_segment(aes(x=filter, xend=filter, y=0, yend=value)) +
  xlab("Filter")+ ylab("% samples with at least 1 candidate gene")+ RD_theme
#ggsave('non_zero_samples_expression_plot.pdf', non_zero_samples_plot, path='/srv/scratch/restricted/rare_diseases/analysis/manuscript/figures/', width=4.5, height=6)


# Splicing non zero individual in %
spli_non_zero_samples=colSums(sample_outlier_filter_withsgl.df[,3:(ncol(sample_outlier_filter_withsgl.df)-3)] != 0)/length(sample_with_data) *100
spli_non_zero_samples=melt(spli_non_zero_samples)
spli_non_zero_samples$filter=c(1:8)



spli_non_zero_samples_plot=ggplot(spli_non_zero_samples, aes(x=as.factor(filter), y=value)) + 
  geom_point(size=3) + 
  geom_segment(aes(x=filter, xend=filter, y=0, yend=value)) +
  xlab("Filter")+ ylab("% samples with at least 1 candidate gene")+ RD_theme
#ggsave('non_zero_samples_plot.pdf', non_zero_samples_plot, path='/srv/scratch/restricted/rare_diseases/analysis/manuscript/figures/', width=4, height=6)



## combined_plot
sup_percent_nonzero_plot=ggdraw()+draw_plot(exp_non_zero_samples_plot, 0,0,1/2,1)+
	draw_plot(spli_non_zero_samples_plot, 1/2,0,1/2,1)+
	draw_plot_label(c('A', 'B'), c(0,1/2), c(1,1), size = 20)

sup_percent_nonzero_plot

pdf('/srv/scratch/restricted/rare_diseases/analysis/manuscript/figures/sup_percent_nonzero_plot_v2.pdf', w=11, h=8)
sup_percent_nonzero_plot
dev.off()

