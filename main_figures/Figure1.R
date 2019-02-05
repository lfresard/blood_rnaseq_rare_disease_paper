#!/bin/R
# LF
# November 2017


# Code to generate figure 1.
# Figure 1A: pie chart of broad disease categories
# Figure 1B: Expression in blood of disease genes
# Figure 1C: %genes expressed in blood when gene is expressed in only one tissue vs more than one tissue
# Figure 1D: expression in blood of genes with pLI>=9


#rm(list=ls())
load(paste0(dir,"/analysis/manuscript/figures_revision/Figure1.RData"))

#--- Libraries

library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(dplyr)

#--- Functions

source('/users/lfresard/repos/rare_disease/scripts/manuscript_analyses/Figures_source.R')


#--- Main

#--- Figures
fsize=15
RD_theme=theme_classic()+
	theme(axis.text.x= element_text(size=fsize),
	axis.text.y= element_text(size=fsize), 
	axis.title = element_text(size = fsize), 
    legend.text = element_text(size = fsize-1), 
    legend.title = element_text(size = fsize), 
    axis.ticks = element_line(size = 0.1))


# Figure 1A
blank_theme <- theme_minimal()+
  theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  plot.title=element_text(size=fsize, face="bold"),axis.text.x= element_text(size=fsize),
	axis.text.y= element_text(size=fsize), 
	axis.title = element_text(size = fsize), 
    legend.text = element_text(size = fsize-1), 
    legend.title = element_text(size = fsize), 
    axis.ticks = element_line(size = 0.1)
  )

# grey scale plot
y.breaks <- cumsum(categ_case$counts) - categ_case$counts/2
categ_case =categ_case %>% mutate(cumulative=cumsum(counts),midpoint=cumulative-(counts/2))
categ_case_greyscale_plot= ggplot(categ_case, aes(x="", y=counts, fill=disease_category))+
	geom_bar(width = 1, stat = "identity")+
	coord_polar("y", start=0)+ 
	scale_fill_grey(start = 0.6, end = 1) +  
	blank_theme + theme(legend.position="top")+
	labs(fill="Disease category")+
	scale_y_continuous(
        breaks=y.breaks,   # where to place the labels
        labels=categ_case$disease_category)+# the labels
    theme(legend.position="", axis.text = element_text(colour="black", size = 10))+ 
    geom_text(aes(x=1.2, y=midpoint, label=counts), color="black", fontface="bold", size=3)


ggsave('Figure1A_nf_greyscale.pdf', categ_case_greyscale_plot, path='./figures_revision/', width=6, height=6)



# Figure 1B
disease_colors=brewer.pal(12,"Paired")

col_panelB=c(grey.colors(n=11,start = 0.6, end = 1)[1:4],disease_colors[10])

pct_in_blood_df=pct_in_blood_df %>% select(PCT_20,PCT_50,PCT_90,DISEASE)
pct_in_blood_df$DISEASE=rownames(pct_in_blood_df)
pct_in_blood_df$DISEASE=rownames(pct_in_blood_df)


disease_barplot3_greyscale= ggplot(pct_in_blood_bin_median_df, aes(x=bin, y=pct_genes, fill=DISEASE))+
	geom_bar(stat="identity", position="dodge", color="black")+ 
	scale_fill_manual(values=col_panelB,name="",guide = guide_legend(nrow=2)) +
	labs(x='Average TPM', y='% Genes expressed in blood')+
	geom_text(stat='identity',aes(label=pct_in_blood_bin_df$DISEASE),vjust=1.6,hjust=1 ,color="white", size=fsize/3, position = position_dodge(width=1))+
	RD_theme+
	coord_flip()+
	#geom_text(aes(x=bin, y=6, label=pct_in_blood_bin_df$DISEASE), color="white", fontface="bold", size=3.3, position = position_dodge())+
	theme(legend.position="")
ggsave('Figure1B_nf_greyscale.pdf', disease_barplot3_greyscale, path='./figures_revision/', width=6, height=6)


disease_barplot3_greyscale_cumsum= ggplot(pct_in_blood_bin_median_df3, aes(x=bin, y=pct_genes, fill=DISEASE))+
	geom_bar(stat="identity", position="dodge", color="black")+ 
	scale_fill_manual(values=col_panelB,name="",guide = guide_legend(nrow=2)) +
	labs(x='Median TPM', y='% Genes expressed in blood')+
	geom_text(stat='identity',aes(label=pct_in_blood_bin_median_df3$DISEASE),vjust=1,hjust=1 ,color="white", size=fsize/3, position = position_dodge(width=1))+
	RD_theme+
	coord_flip()+
	#geom_text(aes(x=bin, y=6, label=pct_in_blood_bin_df$DISEASE), color="white", fontface="bold", size=3.3, position = position_dodge())+
	theme(legend.position="")
disease_barplot3_greyscale_cumsum
ggsave('Figure1B_nf_greyscale_cumsum.pdf', disease_barplot3_greyscale_cumsum, path='./figures_revision/', width=6, height=6)


# Figure 1
combined_plots_greyscale=ggdraw()+draw_plot(categ_case_greyscale_plot, 0,0,1/2,1)+
	draw_plot(disease_barplot3_greyscale_cumsum, 1/2,1/3,1/2,2/3)+
	#draw_plot(exp_multtissues_genes_plot_greyscale, 0,0,1/2,1/2)+
	#draw_plot(exp_blood_pLI_plot, 1/2,0,1/2,1/2)+

	draw_plot_label(c('A', 'B'), c(0,1/2), c(1,1), size = 15)

#combined_plots
pdf('./figures_revision/Figure1.pdf', w=9, h=9)
combined_plots_greyscale
dev.off()


