#!/bin/R

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


# load required functions
source('Figures_source.R')
#source('/users/lfresard/repos/rare_disease/scripts/manuscript_analyses/Figures_source.R')


##--- MAIN

# List junction files 
junction_dir=paste(dir,"/data/splicing/juncfiles/all_unfiltjunc/",sep="")
junc_suffix=".SJ.out_gene_info.tsv"

# Read rare disease samples metadata
metadata = paste(dir,"/data/metadata/2018_12_02_Rare_Disease_Metadata.tsv",sep="")
metadata = read_tsv(metadata) %>% filter(in_freeze=="yes",status=="PASSED", is_RD=="yes",source_of_RNA=="Blood")
RD_samples=metadata$sample_id

# Get metadata info for pivus and dgn
meta_pivus2=paste(dir,"/data/metadata/Pivus_expCovariatesAll_unscaled.txt",sep="")
meta_pivus2=read_tsv(meta_pivus2)%>% filter(Age==0) 
pivus_samples=meta_pivus2$Sample 

dgn_junc=list.files(path=junction_dir, pattern=glob2rx("LD*.SJ.out_gene_info.tsv"))
dgn_samples=str_extract(dgn_junc, "LD[0-9]+")


# get all samples id
samples=c(RD_samples,pivus_samples,dgn_samples)

# select junc files corresponding to samples
junc_files=paste0(junction_dir,samples,junc_suffix)
samples=samples[(!file.size(junc_files) == 0)]
junc_files=junc_files[(!file.size(junc_files) == 0)]

#make a dataframe containing all junction counts for selected samples
junc_info=lapply(samples,process_junc_file,junction_dir,junc_suffix)
merged.data.frame = Reduce(function(...) merge(...,  by="junction",all=T), junc_info)
colnames(merged.data.frame)=c("junction",samples)


# Get annotated junctions
annot_junc=paste(dir,"/data/splicing/juncfiles/gencodev19_intronsstartplus1_proteincoding.genenames_uniqjunc.tsv",sep="")
annot_junc=read_tsv(annot_junc, col_names=FALSE)
colnames(annot_junc)=c("chr", "junc_start", "junc_end", "gene")
annot_junctions=paste(annot_junc$chr, annot_junc$junc_start, annot_junc$junc_end, annot_junc$gene, sep="_")

# get number of annotated junction per gene
annot_junc_num_gene=melt(table(str_extract(annot_junc$gene, "ENSG[0-9]+")))

# Filter merged data frame for junctions that are annotated
annot.junc=merged.data.frame %>% filter(junction %in% annot_junctions)

# filter out junction for which less than 20% of samples with at least 5 reads.
annot.junc[is.na(annot.junc)]<-0 # put NA junctions to 0
annot.junc.filt=annot.junc[rowSums(annot.junc[,2:ncol(annot.junc)]>=5) >(length(samples) *20/100),]

# get genes for which we observe this coverage
nb_junction_cov=melt(table(str_extract(annot.junc.filt$junction, "ENSG[0-9]+")))
junc_genes=nb_junction_cov$Var1
nb_junction_cov=nb_junction_cov %>% left_join(annot_junc_num_gene, by="Var1")
names(nb_junction_cov)=c("gene", "covered", "annotated")
nb_junction_cov=nb_junction_cov%>% mutate(percent_covered=covered/annotated)



# Read in disease gene lists
disease_genes_dir=paste(dir,"/data/candidate_gene_lists/",sep="")
gene_lists=list.files(path=disease_genes_dir, pattern='genelist.txt')

disease_lists=lapply(gene_lists, read_in_disease_lists, path_to_file=disease_genes_dir)
names(disease_lists)=unlist(strsplit(gene_lists, "_genelist.txt"))

# get ensembl names for all gene lists
disease_lists_ens=lapply(disease_lists,check_gene_names)



# get overlap between genes with coverage and disease lists
disease_junc=data.frame(OMIM=length( disease_lists_ens$OMIM[disease_lists_ens$OMIM %in% junc_genes])/length(disease_lists_ens$OMIM), Neurology=length( disease_lists_ens$Neurology[disease_lists_ens$Neurology %in% junc_genes])/length(disease_lists_ens$Neurology),Musculoskeletal=length( disease_lists_ens$Skeletal_disorders[disease_lists_ens$Skeletal_disorders %in% junc_genes])/length(disease_lists_ens$Skeletal_disorders), Hematology=length( disease_lists_ens$Hematology[disease_lists_ens$Hematology %in% junc_genes])/length(disease_lists_ens$Hematology), Ophthalmology=length( disease_lists_ens$Ophtalmology[disease_lists_ens$Ophtalmology %in% junc_genes])/length(disease_lists_ens$Ophtalmology))
disease_junc.m=melt(disease_junc)
disease_junc.m$variable=factor(disease_junc.m$variable, levels=c("Hematology", "Ophthalmology","Musculoskeletal", "Neurology", "OMIM"))



##-------------------------
# Genes expressed in more than 1 tissue and expression in Blood
# We looked at GTEx v7 RPKM for this plot

# Read in % results
exp_multtissues_file="/srv/scratch/restricted/rare_diseases/data/gtex_tissues_exp/Blood_genes_comp_othertissues_exp.txt"
exp_multtissues=read.table(exp_multtissues_file, header=T, sep="\t")


# Make data frame out of results
multis_df=data.frame(tissues=c('1 tissue', '>1 tissue'), 
	overlap_blood=c(exp_multtissues$gene_in_blood[exp_multtissues$tissue_number==1]*100, sum(exp_multtissues$gene_in_blood[2:nrow(exp_multtissues)])*100),
	N_genes=c(exp_multtissues$gene_number[exp_multtissues$tissue_number==1],sum(exp_multtissues$gene_number[2:nrow(exp_multtissues)])))

multis_df$tissues=factor(multis_df$tissues, levels=rev(levels(multis_df$tissues)))


## Look at gene level
exp_multtissues_genes_file="/srv/scratch/restricted/rare_diseases/data/gtex_tissues_exp/Number_of_tissue_expression.txt"
exp_multtissues_genes=as.data.frame(read.table(exp_multtissues_genes_file, header=F, sep="\t"))
colnames(exp_multtissues_genes)=c("ensgene", "tissues")
exac = as.data.frame(read.table('/users/lfresard/NumberOfIndividuals_Impact/enrichment_disease-conservation_database/annotations/EXAC/forweb_cleaned_exac_r03_march16_z_data_pLI.txt', sep = '\t', stringsAsFactors = F, header = T))

# stratify analysis between high/low pLIs
exp_multtissues_genes=exp_multtissues_genes %>% 
	inner_join(grch37, by="ensgene") %>% 
	select(ensgene, tissues, symbol) %>% 
	distinct %>% inner_join(exac, by=c("symbol"="gene")) %>% 
	select(ensgene,tissues,syn_z,mis_z,lof_z,pLI) %>% 
	distinct %>% mutate(multi=ifelse(tissues==1,"no", "yes"))%>%select(-one_of(c("tissues", "pLI")))


exp_multtissues_genes.m=melt(exp_multtissues_genes, id.vars = c("ensgene", "multi"))


## expression of extreme lofz and misz in blood

# get average expression in blood
tpms_df_avg=data.frame(ensgene=tpms_df$ensgene, avg_exp=rowMeans(tpms_df[,2:ncol(tpms_df)]))
tpms_df_avg=tpms_df_avg %>%
	inner_join(grch37, by="ensgene") %>% filter(biotype == "protein_coding")%>%
	select(ensgene, symbol,avg_exp) %>% 
	inner_join(exac, by=c("symbol"="gene")) %>% 
	select(ensgene,avg_exp,syn_z,mis_z,lof_z, pLI) %>% distinct %>% filter(pLI>=.9)

#breaks=c(0, 0.1,1,10,100,10000)
breaks=c(0, 0.1,10,10000,100000)
tpms_df_avg=transform(tpms_df_avg, bin = cut(avg_exp, breaks, right=F))



# Save data
save.image(file = paste0(dir,"/analysis/manuscript/figures_revision/FigureS2.in.RData"))
