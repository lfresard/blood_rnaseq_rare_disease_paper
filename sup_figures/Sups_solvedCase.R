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


# select results for rars2_ samples
rars2_1.df=exp_outlier_number.df %>% filter(sample=="RD182")
rars2_2.df=exp_outlier_number.df %>% filter(sample=="RD194")


# remove RIVER scores
#rars2_1.df=rars2_1.df[-9,]
#rars2_2.df=rars2_2.df[-9,]


# Annotate with filter category
rars2_1.df$category=c("Rare variant", rep("Rare variant and filters",3),"HPO", "Expression outlier", rep("Expression outlier and filters",7))
rars2_1.df$category=factor(rars2_1.df$category, levels=c("Expression outlier and filters","Expression outlier","HPO","Rare variant and filters","Rare variant"))

rars2_2.df$category=c("Rare variant", rep("Rare variant and filters",3),"HPO", "Expression outlier", rep("Expression outlier and filters",7))
rars2_2.df$category=factor(rars2_2.df$category, levels=c("Expression outlier and filters","Expression outlier","HPO","Rare variant and filters","Rare variant"))

# Get genes for each filter

rars2_1.genes=get_genes_filters("RD182")
rars2_2.genes=get_genes_filters("RD194")


# remove sample name 
rars2_1.genes=rars2_1.genes[-1]
rars2_2.genes=rars2_2.genes[-1]

# look if causal gene

rars2_1.df$is_causal_in=as.factor(unlist(lapply(rars2_1.genes, is_causal,"rars2_")))
rars2_2.df$is_causal_in=as.factor(unlist(lapply(rars2_2.genes, is_causal,"rars2_")))


# plot results
rars2_1.dotchart=ggdotchart(rars2_1.df, x = "filter", y = "value",color = "category",shape='is_causal_in',
	palette=c("#D0C6B1","#FFD800", "#587058", "#587498", "#E86850"),
	sorting = "descending",add = "segments",rotate = TRUE,
	group = "category",dot.size =20,label = round(rars2_1.df$value),
	font.label = list(color = "white", size = 18, vjust = 0.5), ggtheme = theme_pubr(),
	add.params = list(color = "category", size = 2), 
	#legend=c(0.8,0.5), ylab="Number of candidate genes", xlab="Filter")#+RD_theme
	legend="", ylab="Number of candidate genes", xlab="Filter")#+RD_theme
rars2_1.dotchart

rars2_2.dotchart=ggdotchart(rars2_2.df, x = "filter", y = "value",color = "category",shape='is_causal_in',
	palette=c("#D0C6B1","#FFD800", "#587058", "#587498", "#E86850"),
	sorting = "descending",add = "segments",rotate = TRUE,
	group = "category",dot.size =20,label = round(rars2_2.df$value),
	font.label = list(color = "white", size = 18, vjust = 0.5), ggtheme = theme_pubr(),
	add.params = list(color = "category", size = 2), 
	#legend=c(0.8,0.5), ylab="Number of candidate genes", xlab="Filter")#+RD_theme
	legend="", ylab="Number of candidate genes", xlab="Filter")#+RD_theme
rars2_2.dotchart


rars2_1_rank.df <- data.frame(filter=c("Expression outliers - Intitial count", "pLI over 0.9","ASE in gene", "RIVER score over 0.85","HPO term matches gene", "Rare variant within 10kb","RV + CADD score over 10","RV + CADD score over 10 + HPO term match"), rank=c(666,NA,NA,NA,7,261,39,10))
rars2_2_rank.df <- data.frame(filter=c("Expression outliers - Intitial count", "pLI over 0.9","ASE in gene", "RIVER score over 0.85","HPO term matches gene", "Rare variant within 10kb","RV + CADD score over 10","RV + CADD score over 10 + HPO term match"), rank=c(130,NA,NA,NA,6,39,7,2))


rars2_1_rank.df$filter=factor(rars2_1_rank.df$filter, levels=c( "RV + CADD score over 10 + HPO term match","RV + CADD score over 10","Rare variant within 10kb","HPO term matches gene", "RIVER score over 0.85","ASE in gene","pLI over 0.9", "Expression outliers - Intitial count"))
rars2_2_rank.df$filter=factor(rars2_2_rank.df$filter, levels=c( "RV + CADD score over 10 + HPO term match","RV + CADD score over 10","Rare variant within 10kb","HPO term matches gene", "RIVER score over 0.85","ASE in gene","pLI over 0.9", "Expression outliers - Intitial count"))



rars2_1_rank.plot=ggplot(rars2_1_rank.df, aes(x=filter, y=rank)) + scale_y_continuous(breaks = pretty(rars2_1_rank.df$rank, n = 666))+#RD_theme+
  geom_point(col="#D0C6B1", size=20) +   # Draw points
  geom_segment(x=1:8, xend=1:8, y=1, yend=666, linetype="dashed", size=0.1) + labs(y="zscore rank", x="") +theme(axis.text.y=element_blank())+ # Draw dashed lines
  #labs(title="Dot Plot", subtitle="Make Vs Avg. Mileage", caption="source: mpg") +  
  coord_flip()

rars2_2_rank.plot=ggplot(rars2_2_rank.df, aes(x=filter, y=rank)) + scale_y_continuous(breaks = pretty(rars2_1_rank.df$rank, n = 130))+#RD_theme+
  geom_point(col="#D0C6B1", size=20) +   # Draw points
  geom_segment(x=1:8, xend=1:8, y=1, yend=130, linetype="dashed", size=0.1) + labs(y="zscore rank", x="") +theme(axis.text.y=element_blank())+ # Draw dashed lines
  #labs(title="Dot Plot", subtitle="Make Vs Avg. Mileage", caption="source: mpg") +  
  coord_flip()

pdf("figures_revision/rars2_1.rank.pdf", height=5, width=15)
rars2_1_rank.plot
dev.off()

pdf("figures_revision/rars2_2.rank.pdf", height=5, width=15)
rars2_2_rank.plot
dev.off()

## combined_plot
rars2_1_big_plot=ggdraw()+draw_plot(rars2_1.dotchart, 0,0,2/3,1)+
	draw_plot(rars2_1_rank.plot, 2/3,0,1/3,1)+
	draw_plot_label(c('A'), c(0), c(1), size = 20)

pdf('figures_revision/rars2_1_big_plot.pdf', w=10, h=12)
rars2_1_big_plot
dev.off()




## combined_plot
rars2_2_big_plot=ggdraw()+draw_plot(rars2_2.dotchart, 0,0,2/3,1)+
	draw_plot(rars2_2_rank.plot, 2/3,0,1/3,1)+
	draw_plot_label(c('A'), c(0), c(1), size = 20)

pdf('figures_revision/rars2_2_big_plot.pdf', w=10, h=12)
rars2_2_big_plot
dev.off()










# select results for rars2_ samples
rars2_1.df=exp_outlier_number.df %>% filter(sample=="RD059")
rars2_2.df=exp_outlier_number.df %>% filter(sample=="RD062")


# remove RIVER scores
#rars2_1.df=rars2_1.df[-9,]
#rars2_2.df=rars2_2.df[-9,]


# Annotate with filter category
rars2_1.df$category=c("Rare variant", rep("Rare variant and filters",3),"HPO", "Expression outlier", rep("Expression outlier and filters",7))
rars2_1.df$category=factor(rars2_1.df$category, levels=c("Expression outlier and filters","Expression outlier","HPO","Rare variant and filters","Rare variant"))

rars2_2.df$category=c("Rare variant", rep("Rare variant and filters",3),"HPO", "Expression outlier", rep("Expression outlier and filters",7))
rars2_2.df$category=factor(rars2_2.df$category, levels=c("Expression outlier and filters","Expression outlier","HPO","Rare variant and filters","Rare variant"))

# Get genes for each filter

rars2_1.genes=get_genes_filters("RD059")
rars2_2.genes=get_genes_filters("RD062")


# remove sample name 
rars2_1.genes=rars2_1.genes[-1]
rars2_2.genes=rars2_2.genes[-1]

# look if causal gene

rars2_1.df$is_causal_in=as.factor(unlist(lapply(rars2_1.genes, is_causal,"RARS2")))
rars2_2.df$is_causal_in=as.factor(unlist(lapply(rars2_2.genes, is_causal,"RARS2")))


# plot results
rars2_1.dotchart=ggdotchart(rars2_1.df, x = "filter", y = "value",color = "category",shape='is_causal_in',
	palette=c("#D0C6B1","#FFD800", "#587058", "#587498", "#E86850"),
	sorting = "descending",add = "segments",rotate = TRUE,
	group = "category",dot.size =20,label = round(rars2_1.df$value),
	font.label = list(color = "white", size = 18, vjust = 0.5), ggtheme = theme_pubr(),
	add.params = list(color = "category", size = 2), 
	#legend=c(0.8,0.5), ylab="Number of candidate genes", xlab="Filter")#+RD_theme
	legend="", ylab="Number of candidate genes", xlab="Filter")#+RD_theme
rars2_1.dotchart

rars2_2.dotchart=ggdotchart(rars2_2.df, x = "filter", y = "value",color = "category",shape='is_causal_in',
	palette=c("#D0C6B1","#FFD800", "#587058", "#587498", "#E86850"),
	sorting = "descending",add = "segments",rotate = TRUE,
	group = "category",dot.size =20,label = round(rars2_2.df$value),
	font.label = list(color = "white", size = 18, vjust = 0.5), ggtheme = theme_pubr(),
	add.params = list(color = "category", size = 2), 
	#legend=c(0.8,0.5), ylab="Number of candidate genes", xlab="Filter")#+RD_theme
	legend="", ylab="Number of candidate genes", xlab="Filter")#+RD_theme
rars2_2.dotchart


rars2_2_rank.df <- data.frame(filter=c("Expression outliers - Intitial count", "pLI over 0.9","ASE in gene", "RIVER score over 0.85","HPO term matches gene", "Rare variant within 10kb","RV + CADD score over 10","RV + CADD score over 10 + HPO term match"), rank=c(331,NA,NA,NA,33,12,7,2))


#rars2_1_rank.df$filter=factor(rars2_1_rank.df$filter, levels=c( "RV + CADD score over 10 + HPO term match","RV + CADD score over 10","Rare variant within 10kb","HPO term matches gene", "RIVER score over 0.85","ASE in gene","pLI over 0.9", "Expression outliers - Intitial count"))
rars2_2_rank.df$filter=factor(rars2_2_rank.df$filter, levels=c( "RV + CADD score over 10 + HPO term match","RV + CADD score over 10","Rare variant within 10kb","HPO term matches gene", "RIVER score over 0.85","ASE in gene","pLI over 0.9", "Expression outliers - Intitial count"))



rars2_2_rank.plot=ggplot(rars2_2_rank.df, aes(x=filter, y=rank)) + scale_y_continuous(breaks = pretty(rars2_1_rank.df$rank, n = 331))+#RD_theme+
  geom_point(col="#D0C6B1", size=20) +   # Draw points
  geom_segment(x=1:8, xend=1:8, y=1, yend=331, linetype="dashed", size=0.1) + labs(y="zscore rank", x="") +theme(axis.text.y=element_blank())+ # Draw dashed lines
  #labs(title="Dot Plot", subtitle="Make Vs Avg. Mileage", caption="source: mpg") +  
  coord_flip()

pdf("figures_revision/rars2_2.rank.pdf", height=5, width=15)
rars2_2_rank.plot
dev.off()

## combined_plot
rars2_2_big_plot=ggdraw()+draw_plot(rars2_2.dotchart, 0,0,2/3,1)+
	draw_plot(rars2_2_rank.plot, 2/3,0,1/3,2/3)+
	draw_plot_label(c('B'), c(0), c(1), size = 20)

pdf('figures_revision/rars2_2_big_plot.pdf', w=10, h=12)
rars2_2_big_plot
dev.off()




## combined_plot
rars2_2_big_plot=ggdraw()+draw_plot(rars2_2.dotchart, 0,0,2/3,1)+
	draw_plot(rars2_2_rank.plot, 2/3,0,1/3,1)+
	draw_plot_label(c('A'), c(0), c(1), size = 20)

pdf('figures_revision/rars2_2_big_plot.pdf', w=10, h=12)
rars2_2_big_plot
dev.off()


pdf('figures_revision/rars2_1.dotchart.pdf', w=10, h=12)
rars2_1.dotchart
dev.off()

# Save data
save.image(file = paste(dir,"/analysis/manuscript/figures_revision/Figure6.in.RData",sep=""))

