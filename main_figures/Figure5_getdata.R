library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(reshape2)
library(data.table)
library(annotables)
library(ggpubr)
library(cowplot)



is_causal=function(item,symbol){
	if (any(item==symbol)==TRUE){
		causal=1
	}else{
		causal=0
	}
	return(causal)
}


# select results for MECR samples
mecr1.df=exp_outlier_number.df %>% filter(sample=="RD182")
mecr2.df=exp_outlier_number.df %>% filter(sample=="RD194")


# remove RIVER scores
#mecr1.df=mecr1.df[-9,]
#mecr2.df=mecr2.df[-9,]


# Annotate with filter category
mecr1.df$category=c("Rare variant", rep("Rare variant and filters",3),"HPO", "Expression outlier", rep("Expression outlier and filters",7))
mecr1.df$category=factor(mecr1.df$category, levels=c("Expression outlier and filters","Expression outlier","HPO","Rare variant and filters","Rare variant"))

mecr2.df$category=c("Rare variant", rep("Rare variant and filters",3),"HPO", "Expression outlier", rep("Expression outlier and filters",7))
mecr2.df$category=factor(mecr1.df$category, levels=c("Expression outlier and filters","Expression outlier","HPO","Rare variant and filters","Rare variant"))

# Get genes for each filter

mecr1.genes=get_genes_filters("RD182")
mecr2.genes=get_genes_filters("RD194")


# remove sample name 
mecr1.genes=mecr1.genes[-1]
mecr2.genes=mecr2.genes[-1]

# look if causal gene

mecr1.df$is_causal_in=as.factor(unlist(lapply(mecr1.genes, is_causal,"MECR")))
mecr2.df$is_causal_in=as.factor(unlist(lapply(mecr2.genes, is_causal,"MECR")))


# plot results
mecr1.dotchart=ggdotchart(mecr1.df, x = "filter", y = "value",color = "category",shape='is_causal_in',
	palette=c("#D0C6B1","#FFD800", "#587058", "#587498", "#E86850"),
	sorting = "descending",add = "segments",rotate = TRUE,
	group = "category",dot.size =20,label = round(mecr1.df$value),
	font.label = list(color = "white", size = 18, vjust = 0.5), ggtheme = theme_pubr(),
	add.params = list(color = "category", size = 2), 
	#legend=c(0.8,0.5), ylab="Number of candidate genes", xlab="Filter")#+RD_theme
	legend="", ylab="Number of candidate genes", xlab="Filter")#+RD_theme
kctd7.dotchart

