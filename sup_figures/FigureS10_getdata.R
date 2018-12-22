#!/bin/R

# look at the number of outliers depending on the junction number

#--- Libraries
library(readr)
library(dplyr)
library(reshape2)
library(cowplot)
library(annotables)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(tximport)
library(stringr)
library(gtools)
library(ggpubr)



# Master directory
dir = Sys.getenv('RARE_DIS_DIR')


# Get annotated junctions
annot_junc=paste(dir,"/data/splicing/juncfiles/gencodev19_intronsstartplus1_proteincoding.genenames_uniqjunc.tsv", sep="")
annot_junc=read_tsv(annot_junc, col_names=FALSE)
colnames(annot_junc)=c("chr", "junc_start","junc_end", "gene")
annot_junctions=paste(annot_junc$chr, annot_junc$junc_start, annot_junc$junc_end, annot_junc$gene, sep="_")

#number of annotated junction per gene
annot_junc_num_gene=melt(table(str_extract(annot_junc$gene, "ENSG[0-9]+")))

# Get splicing outlier
zscores_blood =as.data.frame(read.table(paste(dir,"/analysis/outlier_analysis/splicing/for_freeze/RD_freeze_junc_outliers_PCAratio_blood_nf_PIVUS_withPIVsamples_withmisdata.txt",sep=""), sep="\t", header=T))

# get junction cluster id
ratio_file=paste(dir,"/analysis/outlier_analysis/splicing/for_freeze/RD_PIVUS_freeze_junc_ratios.txt", sep="")
ratios=as.data.frame(read.table(ratio_file,  header=T, sep="\t"))
ratios=ratios %>% select(chr, junction_start, junction_end, gene, annotation_status, Cluster)

# Read rare disease samples metadata
metadata = paste(dir,"/data/metadata/2018_06_12_Rare_Disease_Metadata.tsv",sep="")
metadata = read_tsv(metadata) %>% filter(in_freeze=="yes") %>% filter(source_of_RNA=="Blood")
cases=metadata %>% filter(affected_status=="Case") %>% select(sample_id) %>% pull

# Combine outlier info with annotated data
z2=zscores_blood %>% mutate(affected_status=ifelse(variable %in% cases, "Case", "Control"))%>%
	left_join(ratios, by=c("chr"="chr", "junc_start"="junction_start", "junc_end"="junction_end", "gene"="gene", "annotation_status"="annotation_status")) %>% 
	filter(absZ>=2)  %>% 
	select(variable,gene,affected_status, Cluster)%>%
	filter(!duplicated(Cluster)) %>%
	mutate(gene=str_extract(gene,"ENSG[0-9]+")) %>%
	group_by(variable, gene,affected_status) %>%
	summarize(outliers=n()) %>%
	ungroup %>%
	left_join(annot_junc_num_gene, by=c("gene"="Var1")) %>% 
	group_by(gene,affected_status, value) %>%
	summarize(average=mean(outliers))  %>% mutate(juncbin=cut(value, breaks=c(0,  5,20, 50,  153), right=TRUE))

save.image(file = paste(dir,"/data/FigureS10.in.RData",sep=""))
