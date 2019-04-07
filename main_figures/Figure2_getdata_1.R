#Libraries
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(sva)
library(cowplot)
library(purrr)
library(broom)
library(reshape2)
library(ggfortify)
library(gridExtra)
library(data.table)
library(annotables)
library(corrplot)
library(RColorBrewer)
library(qvalue) 
library(ggpubr)
library(rlist)
library(ggbeeswarm)
library(foreach)
require(doMC)

rm(list=ls())
dir = Sys.getenv('RARE_DIS_DIR')


load(file=paste0(dir,"/analysis/manuscript/figures_revision/Figure2.in.RData"))

source('/users/lfresard/repos/rare_disease/scripts/manuscript_analyses/Figures_source.R')

# Get metadata
metadata = "/srv/scratch/restricted/rare_diseases/data/metadata/2018_12_18_Rare_Disease_Metadata.tsv"
metadata = read_tsv(metadata) %>% filter(in_freeze=="yes") %>% filter(source_of_RNA=="Blood",sequencing_status=="PASSED")
metadata = read_tsv(metadata) %>% filter(sequencing_status=="PASSED", is_RD=="yes",source_of_RNA=="Blood")
sample_with_data=metadata %>% filter(variant_data=="exome" | variant_data=="genome")%>% select(sample_id) %>% pull


# read in exac data
exac = as.data.frame(read.table('/users/lfresard/NumberOfIndividuals_Impact/enrichment_disease-conservation_database/annotations/EXAC/forweb_cleaned_exac_r03_march16_z_data_pLI.txt', sep = '\t', stringsAsFactors = F, header = T))
exac=exac %>% inner_join(grch37, by=c("gene"="symbol")) %>% select(ensgene,syn_z, mis_z, lof_z, pLI) 


# get gene strand info
gene_strand=grch37 %>% select(ensgene,strand)



# Read in expression data

exp_outlier="/srv/scratch/restricted/rare_diseases/analysis/outlier_analysis/expression_level/outliers_zscore_pair_spline_07dec18.txt"
exp_outlier=read_tsv(exp_outlier)
# Read in expression outlier with RV
RV_outlier_exp_10kb="/srv/scratch/restricted/rare_diseases/analysis/outlier_analysis/expression_level/outliers_rare_var_combined_10kb_07dec18.txt"
RV_outlier_exp_10kb=read_tsv(RV_outlier_exp_10kb, col_names=FALSE)
colnames(RV_outlier_exp_10kb)=c("sample_id", "gene", "expressionZScore", "variantPos", "gnomadMAF", "raw_cadd", "phred_cadd", "variant_gene_pos")

RV_outlier_exp_10kb=RV_outlier_exp_10kb %>% left_join(grch37, by=c("gene"="ensgene"), all.x=T)%>%select(sample_id,gene,expressionZScore,gnomadMAF,phred_cadd,variant_gene_pos,variantPos, start,end)
RV_outlier_exp_10kb=RV_outlier_exp_10kb%>% mutate(distance=ifelse((variantPos<=start&variantPos<=end)|(variantPos>=start&variantPos>=end),start-variantPos,0))

# make a column with singleton included or expluded form the analysis (window method)
RV_outlier_exp_10kb =RV_outlier_exp_10kb %>% mutate(gnomAD_AF_withsgl = replace(gnomadMAF, gnomadMAF==".", 0))

# format variant data so that only the max observed allele frequency is kept for variant where several frequencies are available (window method)
RV_outlier_exp_10kb$gnomAD_AF_filt_withsgl=sapply(RV_outlier_exp_10kb$gnomAD_AF_withsgl,process_variant)

RV_outlier_exp_10kb=RV_outlier_exp_10kb%>% left_join(gene_strand, by=c("gene"="ensgene"))
RV_outlier_exp_10kb=RV_outlier_exp_10kb%>% mutate(variant_gene_pos2=case_when(variant_gene_pos=="downstream" & strand==1 ~ "downstream", variant_gene_pos=="downstream" &  strand==-1 ~ "upstream", variant_gene_pos=="upstream" & strand==1 ~ "upstream", variant_gene_pos=="upstream" & strand==-1 ~"downstream"))



# filter for variant of interest
rv_outlier_annot=rv_annotation %>% left_join(RV_outlier_exp_filt_withsgl_10kb, by=c("gene", "variantPos", "indiv_id"="sample_id"))



#filter for RV with MAF<=0.1% (window method)
maf=0.001
RV_outlier_exp_filt_withsgl_10kb=RV_outlier_exp_10kb %>% filter(!is.na(gnomAD_AF_filt_withsgl)) %>% filter(gnomAD_AF_filt_withsgl <=maf)

# Read in number of RV in genes
RV_10kb="/srv/scratch/restricted/rare_diseases/analysis/outlier_analysis/expression_level/rare_var_protein_coding_genes_10kb_07dec18.txt"
RV_10kb=read_tsv(RV_10kb, sep="\t", col_names=F)
colnames(RV_10kb)=c("sample_id", "gene",  "variantPos", "gnomadMAF", "raw_cadd", "phred_cadd","variant_gene_pos")


# make a column with singleton included or expluded form the analysis (window method)
RV_10kb =RV_10kb %>% mutate(gnomAD_AF_withsgl = replace(gnomadMAF, gnomadMAF==".", 0))


# format variant data so that only the max observed allele frequency is kept for variant were several frequencies are available (window method)
RV_10kb$gnomAD_AF_filt_withsgl=sapply(RV_10kb$gnomAD_AF_withsgl,process_variant)


#filter for RV with MAF<=0.1% (window method)
maf=0.001
	# Including singletons
RV_withsgl_10kb=RV_10kb %>% filter(!is.na(gnomAD_AF_filt_withsgl)) %>% filter(gnomAD_AF_filt_withsgl <=maf)
RV_withsgl_10kb=RV_withsgl_10kb %>% mutate(gene=str_extract(gene, "ENSG[0-9]+")) %>% left_join(gene_strand, by=c("gene"="ensgene"))

RV_withsgl_10kb=RV_withsgl_10kb %>% mutate(variant_gene_pos2=case_when(variant_gene_pos=="downstream" & strand==1 ~ "downstream", variant_gene_pos=="downstream" &  strand==-1 ~ "upstream", variant_gene_pos=="upstream" & strand==1 ~ "upstream", variant_gene_pos=="upstream" & strand==-1 ~"downstream"))

# filter for variant of interest
RV_withsgl_10kb_annot=rv_annotation %>% left_join(RV_withsgl_10kb, by=c("gene", "variantPos", "indiv_id"="sample_id"))
write.table(RV_withsgl_10kb_annot, "RV_withCADD.txt", quote=F, sep="\t", row.names=T)

# generate HPO terms list for each sample
gene_to_pheno="/srv/scratch/restricted/rare_diseases/data/metadata/18_10_23_ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt"
gene_to_pheno = read_tsv(gene_to_pheno, comment= "#",col_names = FALSE)
colnames(gene_to_pheno)=c("entrez", "symbol", "DiseaseId", "HPO_ID")

pheno_annot="/srv/scratch/restricted/rare_diseases/data/metadata/18_10_23_ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt"
pheno_annot=read_tsv(pheno_annot, comment= "#",col_names = FALSE)
colnames(pheno_annot)=c("HPO_ID", "HPO_name", "entrez", "symbol")

# Read in and process alternative HPO files (generated by parse_hpoterms.py)
parent_terms="/srv/scratch/restricted/rare_diseases/data/metadata/parent_terms_hpo_2018_10_23.txt"
child_terms="/srv/scratch/restricted/rare_diseases/data/metadata/child_terms_hpo_2018_10_23.txt"
altid_terms="/srv/scratch/restricted/rare_diseases/data/metadata/altid_terms_hpo_2018_10_23.txt"

altid_terms=as.data.frame(read.table(altid_terms, header=F, sep="\t"))
colnames(altid_terms)=c("HPO_ID", "ALT_ID")
child_terms=as.data.frame(read.table(child_terms, header=F, sep="\t"))
colnames(child_terms)=c("HPO_ID", "ALT_ID")
parent_terms=as.data.frame(read.table(parent_terms, header=F, sep="\t"))
colnames(parent_terms)=c("HPO_ID", "ALT_ID")


# generates list of HPO alternative IDs from the input files
alt_hpo_ids_list=lapply(altid_terms$HPO_ID,get_hpo_alt_id, altid_terms=altid_terms)
names(alt_hpo_ids_list)=altid_terms$HPO_ID
child_hpo_ids_list=lapply(child_terms$HPO_ID,get_hpo_alt_id, altid_terms=child_terms)
names(child_hpo_ids_list)=child_terms$HPO_ID
parent_hpo_ids_list=lapply(parent_terms$HPO_ID,get_hpo_alt_id, altid_terms=parent_terms)
names(parent_hpo_ids_list)=parent_terms$HPO_ID

# For every ID associated with the disease, get list of alternative Ids, child ids and parent ids
# make a list for each case with HPO terms

cases=metadata %>% filter(affected_status=="Case") %>% select(sample_id) %>% pull

# Get list of HPO terms from metadata file
HPOs=lapply(cases,get_list_hpo, metadata)
names(HPOs)=cases

# select cases with genetic data
cases_withGen=metadata %>% filter(variant_data !="none", affected_status=="Case")%>% select(sample_id) %>% pull
cases_withGen=cases_withGen[!cases_withGen=="RD008"]

controls_withGen=metadata %>% filter(variant_data !="none", affected_status=="Control")%>% select(sample_id) %>% pull
# Generate extended list of HPO terms for each sample
HPO_extended=lapply(cases, process_sample_ids_hpo, HPO_list=HPOs)
names(HPO_extended)=cases

# generates filtered splicing results without RV
candidates_no_gen_ext=lapply(cases,get_hpo_comp_nogen,zscores_blood,metadata, gene_to_pheno,pheno_annot,HPO_extended)
names(candidates_no_gen_ext)=cases

candidates_zscore_hpo_ext=do.call("rbind", candidates_no_gen_ext)
candidates_zscore_hpo_ext=candidates_zscore_hpo_ext %>% arrange(sample)%>% left_join(metadata, by=c("sample"="sample_id")) %>% select( symbol,overlap_hpo,,max_z ,sample, institution_id, institution)


# generate results filtered splicing with RV
#potential_can=lapply(cases_withGen,get_hpo_comp,RV_outlier_filt_withsgl,metadata, gene_to_pheno,pheno_annot, 10,HPO_extended)
potential_can=lapply(cases_withGen,get_hpo_comp_window,RV_outlier_filt_withsgl,metadata, gene_to_pheno,pheno_annot, 10,HPO_extended)
names(potential_can)=cases_withGen

candidates_zscore_hpo_RV=do.call("rbind", potential_can)
candidates_zscore_hpo_RV=candidates_zscore_hpo_RV %>% arrange(sample)%>% left_join(metadata, by=c("sample"="sample_id")) %>% select( symbol,overlap_hpo,,max_z ,sample, institution_id, institution)


# Read in RIVER scores
river_file="/mnt/lab_data/montgomery/xli6/udn/river_score/UDN_RIVERscore_ByGene.txt"
river=as.data.frame(read.table(river_file, sep="\t", header=T))
	#get list of genes with RIVER score >=0.85 for each sample
gene_river_list=sapply(sample_with_data,get_genes_with_highriver,river)



# Read in ASE data
ase_file="/srv/scratch/restricted/rare_diseases/data/ase/June13/final_merged/ensembleids.ase.filtered.txt"
ase=read_tsv(ase_file, col_names=TRUE)
colnames(ase)=c("sample_id", colnames(ase)[2:length(colnames(ase))])

	# transform ase data 
ase=ase %>% mutate(ensgene=str_extract(Gene, "ENSG[0-9]+")) %>% select(sample_id,Chr,POS, "RefRatio" , "TotalReads","ensgene")%>% filter (! duplicated(ensgene))%>% mutate(RefRatio_trans=ifelse(RefRatio<=0.5,1-RefRatio,RefRatio))%>% filter (! RefRatio_trans==1)
ase=ase%>% group_by(ensgene)%>% mutate(max_RefRatio_trans=max(RefRatio_trans,na.rm=T), median_RefRatio_trans=median(RefRatio_trans, na.rm=T))%>% ungroup
	# get list of genes with ASE for each sample
gene_ase_list=sapply(sample_with_data,get_genes_with_ase,ase)

affected_status_df=metadata%>% select(sample_id,affected_status)

colnames(affected_status_df)=c('sample', 'status')

## Influence of filtering on expression outlier results
zscore_threshold=2
distance_threshold=10000
cadd=10
pli=0.9



sample_exp_outlier_filter_withsgl_10kb_up.df=data.frame(sample_id=cases_withGen,
	"raw"=sapply(cases_withGen, get_outlier_genes_exp, threshold=zscore_threshold, outlier_df.m=exp_outlier, pLI_thre=0, direction="under"),
	"raw.Lof"=sapply(cases_withGen, get_outlier_genes_exp, threshold=zscore_threshold, outlier_df.m=exp_outlier, pLI_thre=pli, direction="under"),
	"raw.HPO"=sapply(cases_withGen,get_hpo_match_nb_exp,outlier_file=exp_outlier,HPOs_list=HPO_extended, pLI_thre=0, direction="under"),
	"raw.ASE"=sapply(cases_withGen, get_outlier_ase_genes, threshold=zscore_threshold, outlier_df=exp_outlier, gene_ase_list=gene_ase_list,pLI_thre=0, direction="under"),
	"raw.river"=sapply(cases_withGen, get_outlier_river_genes, threshold=zscore_threshold, outlier_df=exp_outlier, gene_river_list=gene_river_list,pLI_thre=0, direction="under"),
	"RV"=sapply(cases_withGen, get_outlier_genes_RV_exp_10kb,  outlier_df.m=RV_outlier_exp_filt_withsgl_10kb, CADD=0, pLI_thre=0, direction="under", myvariant_gene_pos="upstream"),
	"RV.LoF.intolerant"=sapply(cases_withGen, get_outlier_genes_RV_exp_10kb,  outlier_df.m=RV_outlier_exp_filt_withsgl_10kb, CADD=0, pLI_thre=pli, direction="under", myvariant_gene_pos="upstream"),
	"RV.CADD"=sapply(cases_withGen, get_outlier_genes_RV_exp_10kb,  outlier_df.m=RV_outlier_exp_filt_withsgl_10kb, CADD=cadd, pLI_thre=0, direction="under", myvariant_gene_pos="upstream"),
	"RV.ASE"=sapply(cases_withGen, get_outlier_genes_RV_exp_10kb_ase,  outlier_df.m=RV_outlier_exp_filt_withsgl_10kb,gene_ase_list=gene_ase_list, CADD=cadd, pLI_thre=0, direction="under", myvariant_gene_pos="upstream"),
	"RV.RIVER"=sapply(cases_withGen, get_outlier_genes_RV_exp_10kb_ase,  outlier_df.m=RV_outlier_exp_filt_withsgl_10kb,gene_ase_list=gene_river_list, CADD=cadd, pLI_thre=0, direction="under", myvariant_gene_pos="upstream"),
	"RV.HPO"=sapply(cases_withGen, get_hpo_match_nb_RV_exp_10kb,  outlier_RV=RV_outlier_exp_filt_withsgl_10kb, HPOs_list=HPO_extended,cadd_thres=0, pLI_thre=0, direction="under", myvariant_gene_pos="upstream"),
	"RV.HPO.CADD"=sapply(cases_withGen, get_hpo_match_nb_RV_exp_10kb,  outlier_RV=RV_outlier_exp_filt_withsgl_10kb, HPOs_list=HPO_extended,cadd_thres=10, pLI_thre=0, direction="under", myvariant_gene_pos="upstream"),
	"RV.HPO.CADD.pLI"=sapply(cases_withGen, get_hpo_match_nb_RV_exp_10kb,  outlier_RV=RV_outlier_exp_filt_withsgl_10kb,HPOs_list=HPO_extended, cadd_thres=cadd, pLI_thre=pli, direction="under", myvariant_gene_pos="upstream"),
	"RV.HPO.CADD.ASE"=sapply(cases_withGen, get_hpo_match_nb_RV_exp_10kb_ase,  outlier_RV=RV_outlier_exp_filt_withsgl_10kb,gene_ase_list=gene_ase_list,HPOs_list=HPO_extended, cadd_thres=cadd, pLI_thre=pli, direction="under", myvariant_gene_pos="upstream"),

	institution=sapply(cases_withGen,get_institution, metadata=metadata),
	affected_status=sapply(cases_withGen,get_affected_status, affected_status_df=affected_status_df),
	technology=sapply(cases_withGen,get_technology, metadata=metadata))

# make a plot with proportions of candidates 

sample_exp_outlier_filter_withsgl_10kb_up.df= sample_exp_outlier_filter_withsgl_10kb_up.df %>% filter(raw!=0)
sample_exp_outlier_filter_withsgl_10kb_up.df$prop_pLI=sample_exp_outlier_filter_withsgl_10kb_up.df$raw.Lof/sample_exp_outlier_filter_withsgl_10kb_up.df$raw
sample_exp_outlier_filter_withsgl_10kb_up.df$prop_HPO=sample_exp_outlier_filter_withsgl_10kb_up.df$raw.HPO/sample_exp_outlier_filter_withsgl_10kb_up.df$raw
sample_exp_outlier_filter_withsgl_10kb_up.df$prop_RIVER=sample_exp_outlier_filter_withsgl_10kb_up.df$raw.river/sample_exp_outlier_filter_withsgl_10kb_up.df$raw
sample_exp_outlier_filter_withsgl_10kb_up.df$prop_ASE=sample_exp_outlier_filter_withsgl_10kb_up.df$raw.ASE/sample_exp_outlier_filter_withsgl_10kb_up.df$raw
sample_exp_outlier_filter_withsgl_10kb_up.df$prop_RV=sample_exp_outlier_filter_withsgl_10kb_up.df$RV/sample_exp_outlier_filter_withsgl_10kb_up.df$raw
sample_exp_outlier_filter_withsgl_10kb_up.df$prop_RV.pLI=sample_exp_outlier_filter_withsgl_10kb_up.df$RV.LoF.intolerant/sample_exp_outlier_filter_withsgl_10kb_up.df$raw
sample_exp_outlier_filter_withsgl_10kb_up.df$prop_CADD=sample_exp_outlier_filter_withsgl_10kb_up.df$RV.CADD/sample_exp_outlier_filter_withsgl_10kb_up.df$raw
sample_exp_outlier_filter_withsgl_10kb_up.df$prop_RV.HPO=sample_exp_outlier_filter_withsgl_10kb_up.df$RV.HPO/sample_exp_outlier_filter_withsgl_10kb_up.df$raw
sample_exp_outlier_filter_withsgl_10kb_up.df$prop_RV.RIVER=sample_exp_outlier_filter_withsgl_10kb_up.df$RV.RIVER/sample_exp_outlier_filter_withsgl_10kb_up.df$raw
sample_exp_outlier_filter_withsgl_10kb_up.df$prop_RV.ASE=sample_exp_outlier_filter_withsgl_10kb_up.df$RV.ASE/sample_exp_outlier_filter_withsgl_10kb_up.df$raw
sample_exp_outlier_filter_withsgl_10kb_up.df$prop_RV.HPO.CADD=sample_exp_outlier_filter_withsgl_10kb_up.df$RV.HPO.CADD/sample_exp_outlier_filter_withsgl_10kb_up.df$raw
sample_exp_outlier_filter_withsgl_10kb_up.df$prop_RV.HPO.CADD.pLI=sample_exp_outlier_filter_withsgl_10kb_up.df$RV.HPO.CADD.pLI/sample_exp_outlier_filter_withsgl_10kb_up.df$raw
sample_exp_outlier_filter_withsgl_10kb_up.df$prop_RV.HPO.CADD.ASE=sample_exp_outlier_filter_withsgl_10kb_up.df$RV.HPO.CADD.ASE/sample_exp_outlier_filter_withsgl_10kb_up.df$raw

sample_exp_outlier_filter_withsgl_10kb_up.df=sample_exp_outlier_filter_withsgl_10kb_up.df%>% select(sample_id,institution,affected_status, technology, prop_pLI, prop_RIVER, prop_ASE, prop_HPO,prop_RV,prop_RV.pLI,prop_CADD,prop_RV.RIVER, prop_RV.HPO, prop_RV.ASE,prop_RV.HPO.CADD,prop_RV.HPO.CADD.pLI,prop_RV.HPO.CADD.ASE)

#melt DF
sample_exp_outlier_filter_withsgl_10kb_up.df.m=melt(sample_exp_outlier_filter_withsgl_10kb_up.df,id.vars=c("sample_id","affected_status", "institution", "technology"))

# select for only cases
sample_exp_outlier_filter_withsgl_10kb_up.df.m.cases=sample_exp_outlier_filter_withsgl_10kb_up.df.m %>% filter(affected_status=="Case")

# apply filters on all genes

exp_outlier_number=sapply(cases_withGen,get_count_filters)
exp_outlier_number=as.data.frame(exp_outlier_number)
exp_outlier_number$filter=rownames(exp_outlier_number)
#exp_outlier_number=exp_outlier_number[-1,]
exp_outlier_number.m=exp_outlier_number %>% gather(variable, value, -filter)

exp_outlier_number.df=data.frame(sample=as.factor(exp_outlier_number.m$variable), value=as.numeric(exp_outlier_number.m$value), filter=unlist(exp_outlier_number.m$filter))
toplot=exp_outlier_number[6:nrow(exp_outlier_number),1:ncol(exp_outlier_number)-1]
toplot=data.matrix(toplot)
percent=t(t(toplot[1:8,]) /toplot[1,]*100)
percent=percent[-1,]
percent.m =melt(percent)

Fig2c=ggplot(percent.m,aes(x=Var1, y=value, fill=Var1))+
	geom_boxplot(color="black", notch=F, show.legend = FALSE)+
	geom_point(size = 0, stroke = 0)+
	scale_fill_manual(values=rep("lightgrey", 7), breaks=c("EXP_OUTLIER_PLI","EXP_OUTLIER_ASE","EXP_OUTLIER_RIVER", "EXP_OUTLIER_HPO","EXP_OUTLIER_RV", "EXP_OUTLIER_RV_CADD", "EXP_OUTLIER_RV_CADD_HPO"),
labels=c(expression("pLI " >="0.9"), "ASE", expression("RIVER " >= "0.85"), "HPO match", expression("Rare variant within 10kb","5 + CADD score  " >= "10"),"4 + 6"), name="Filter")+
	scale_x_discrete(breaks=c("EXP_OUTLIER_PLI","EXP_OUTLIER_ASE","EXP_OUTLIER_RIVER", "EXP_OUTLIER_HPO","EXP_OUTLIER_RV", "EXP_OUTLIER_RV_CADD", "EXP_OUTLIER_RV_CADD_HPO"),
labels=c("1", "2", "3","4","5", "6", "7"))+
	labs(x="Filter", y="Proportion of outlier genes") + RD_theme +guides(fill = guide_legend(override.aes = list(size = 4,shape =c(49,50,51,52,53,54,55))))+theme(legend.position = c(0.8,0.7))
ggsave('fig_2C_prop_candidates_filters_10kb_upstream.pdf', Fig2c, path='/srv/scratch/restricted/rare_diseases/analysis/manuscript/figures_revision/', width=6, height=6)

# plot results
filter_withsglt.exp.case.10kb.plot=ggplot(sample_exp_outlier_filter_withsgl_10kb_up.df.m.cases, aes(x=variable, y=value, fill=variable) )+ 
	geom_boxplot(color="black", notch=F, show.legend = FALSE)+
	geom_point(size = 0, stroke = 0)+
	scale_fill_manual(values=rep("lightgrey", 13), breaks=c("prop_pLI","prop_RIVER","prop_ASE", "prop_HPO","prop_RV", "prop_RV.pLI", "prop_CADD","prop_RV.RIVER","prop_RV.HPO", "prop_RV.ASE", "prop_RV.HPO.CADD","prop_RV.HPO.CADD.pLI","prop_RV.HPO.CADD.ASE"),
		labels=c(expression("pLI " >="0.9"), expression("RIVER " >= "0.85"), "ASE", "HPO match", expression("Rare variant within 10kb","1 + 5","5 + CADD score  " >= "10"), "2 + 5", "4 + 5", "3 + 5", "4 + 7", "1 + 2 + 5 + 7", "2 + 4 + 5 + 7"), name="Filter")+
	scale_x_discrete(breaks=c("prop_pLI","prop_RIVER","prop_ASE", "prop_HPO","prop_RV", "prop_RV.pLI", "prop_CADD","prop_RV.RIVER","prop_RV.HPO", "prop_RV.ASE", "prop_RV.HPO.CADD","prop_RV.HPO.CADD.pLI","prop_RV.HPO.CADD.ASE"),
		labels=c("1", "2", "3","4","5", "6", "7", "8", "9", "10", "11", "12", "13"))+
	labs(x="Filter", y="Proportion of outlier genes") + RD_theme +guides(fill = guide_legend(override.aes = list(size = 4,shape =c(49,50,51,52,53,54,55,56,57,58,59,60,61))))+theme(legend.position = c(0.8,0.7))
filter_withsglt.exp.case.10kb.plot



#exp_outlier_number.df=data.frame(sample=as.factor(exp_outlier_number.m$variable), value=unlist(exp_outlier_number.m$value), filter=unlist(exp_outlier_number.m$filter))
# order factors
exp_outlier_number.df$filter=factor(exp_outlier_number.df$filter, levels=rev(c("RV_10KB", "RV_10KB_CADD", "RV_10KB_HPO", "RV_10KB_CADD_HPO","HPO", "EXP_OUTLIER","EXP_OUTLIER_PLI", "EXP_OUTLIER_HPO","EXP_OUTLIER_ASE","EXP_OUTLIER_RIVER","EXP_OUTLIER_RV","EXP_OUTLIER_RV_CADD","EXP_OUTLIER_RV_CADD_HPO")))

exp_outlier_number_controls=sapply(controls_withGen,get_count_filters)

results_exp_out=mclapply(cases_withGen, get_genes_filters, mc.cores=15)
#results_exp_out_sh=mclapply(cases_withGen, get_genes_filters, mc.cores=5)
names(results_exp_out)=cases_withGen
indiv_id=data.frame(STAN_ID=cases_withGen)
indiv_id=indiv_id %>% left_join(metadata, by=c("STAN_ID"="sample_id"))  %>% select(STAN_ID,institution_id)

# write down all splicing results in serated files for each sample, each sheet being one filter --> does not work
for (i in c(1:length(cases_withGen))){
	print(i)
	for (j in c(2:length(names(results_exp_out[[i]])))) {
		write.table(results_exp_out[[i]][[j]], file=paste0(dir,"/analysis/outlier_analysis/expression_level/",indiv_id[i,2],"_",names(results_exp_out[[i]][j]),"_expression_filtered_results.tsv"), sep="\t", quote=F, row.names=F)
	}
}



# distribution of rare variants around candidates

rv_annotation="/users/xli6/projects/udn/feature_file/udn_variant_annotation.txt"
rv_annotation=read_tsv(rv_annotation)
rv_annotation=rv_annotation %>% mutate(gene=str_extract(gene_id, "ENSG[0-9]+"), variantPos=pos-1)
#test=RV_outlier_exp_filt_withsgl_10kb %>% filter(sample_id=="RD194", gene %in% results_exp_out$RD194$EXP_OUTLIER_RV_CADD_HP$ensgene, variant_gene_pos=="upstream", as.numeric(phred_cadd)>10,gnomAD_AF_filt_withsgl<=0.001)

# 
RV_outlier_exp_filt_withsgl_10kb_upstream_cadd_rv=RV_outlier_exp_filt_withsgl_10kb %>% filter( variant_gene_pos2=="upstream", as.numeric(phred_cadd)>10,gnomAD_AF_filt_withsgl<=0.001)

# Get list of HPO terms from metadata file
CAND=lapply(metadata_udn_cand$sample_id,get_list_candidates, metadata_udn_cand)
names(CAND)=metadata_udn_cand$sample_id

# remove siblings
CAND=list.remove(CAND, c("RD059", "RD182","RD184"))

#remove solved cases
CAND_no_solved=list.remove(CAND, c("RD047", "RD064","RD194", "RD062", "RD159"))

#get genes left after most stringent filter
get_genes_shuffle_real=function(sample, CAND_list){
	OUT_RV_SAMPLE=rv_outlier_annot[rv_outlier_annot$indiv_id==sample & rv_outlier_annot$expressionZScore <=-2 & rv_outlier_annot$variant_gene_pos2=="upstream",] #underexpression only
	exp_out_RV_CADD_HPO=distinct(merge(OUT_RV_SAMPLE,grch37, by.x="gene", by.y="ensgene", all.x=T)[,c("gene", "symbol", "entrez", "expressionZScore")])
	exp_out_RV_CADD_HPO=merge(exp_out_RV_CADD_HPO,pheno_annot, by=c("entrez", "symbol"), all.x=T)
	exp_out_RV_CADD_HPO$hpo_in_input=ifelse(exp_out_RV_CADD_HPO$HPO_ID %in% HPO_extended[[sample]], 1, 0) # annotate rows for which HPO match
	exp_out_RV_CADD_HPO=unique(exp_out_RV_CADD_HPO[,c("gene", "symbol",  "hpo_in_input")])
	exp_out_RV_CADD_HPO=aggregate(exp_out_RV_CADD_HPO$hpo_in_input, by=list(exp_out_RV_CADD_HPO$gene), FUN=sum)
	colnames(exp_out_RV_CADD_HPO)=c("gene", "hpo_match")
	exp_out_RV_CADD_HPO=data.frame(ensgene=unique(exp_out_RV_CADD_HPO$gene[exp_out_RV_CADD_HPO$hpo_match>0]))
	exp_out_RV_CADD_HPO=distinct(merge(exp_out_RV_CADD_HPO, grch37,by="ensgene", all.x=T)[c("ensgene", "symbol")])
	exp_out_RV_CADD_HPO=distinct(merge(exp_out_RV_CADD_HPO, OUT_RV_SAMPLE,by.x="ensgene",by.y="gene", all.x=T)[c("ensgene", "symbol", "expressionZScore")])
	exp_out_RV_CADD_HPO=exp_out_RV_CADD_HPO$symbol 

	if (any(CAND_list[[sample]] %in% exp_out_RV_CADD_HPO)==TRUE){
		candin=1
	}else{candin=0}
	return(candin)
}


results_real=sapply(metadata_udn_cand$sample_id,get_genes_shuffle_real)
results_real=sapply(names(CAND),get_genes_shuffle_real)

unsolved_real=results_real[!names(results_real) %in% c("RD059", "RD062", "RD064", "RD194","RD047", "RD159", "RD062")]


# additional canidtates=sum(unsolved_real)/length(results_real)
registerDoMC(cores = 2)
# create dummy expression zscore datafrmae by shuffling expression value across genes
cand_genes_permut=foreach(i = 1:10000, .combine = rbind) %dopar% {
	
	#zscore_dummy=exp_outlier
	#zscore_dummy$zscore=ave(exp_outlier$zscore,exp_outlier$gene, FUN = sample)

	# select only samples with gendata
	#zscore_dummy=zscore_dummy %>% filter(sample_id %in% metadata_udn_cand$sample_id)
	# merge zscore dummy with RV file

	#RV_dummy=RV_outlier_exp_filt_withsgl_10kb_noZ %>% left_join(zscore_dummy, by=c("gene", "indiv_id"="sample_id"), all.x=T)
	candin=sapply(names(CAND),is_can_in_list)
	sum(candin)/length(names(CAND))*100
}

is_can_in_list=function(sample){
	CAND_list=CAND
	#	 randomize names in CAND list
	names(CAND_list)=sample(names(CAND),replace = FALSE)
	if (any(CAND_list[[sample]] %in% results_exp_out[[sample]][[13]][,2])==TRUE){
		candin=1
	}else{candin=0}
	return(candin)

}
is_can_in_list_real=function(sample){
	if (any(CAND[[sample]]%in%results_exp_out[[sample]][[13]][,2])==TRUE){
		candin=1
	}else{candin=0}
	return(candin)

}
candin=sapply(names(CAND),is_can_in_list_real)
real=sum(candin)/length(names(CAND))*100



candidates_percent=rbind(melt(data.frame(permut=cand_genes_permut)),melt(data.frame(real=real)))

ggplot(candidates_percent, aes(x=variable, y=value, color=variable))+ geom_boxplot()
ggplot(candidates_percent, aes(x=variable, y=value, color=variable))+ geom_jitter()
ggplot(candidates_percent, aes(x=variable, y=value, color=variable)) + geom_quasirandom(varwidth = TRUE)
candidates_percent=data.frame(permut=cand_genes_permut)

is_can_in_list_noexp=function(sample){
	#CAND_list=CAND
	#	 randomize names in CAND list
	#names(CAND_list)=sample(names(CAND),replace = FALSE)
	random_set=sample(results_exp_out[[sample]][[5]][,2], length(results_exp_out[[sample]][[13]][,2]), replace=F)
	if (any(CAND[[sample]] %in% random_set)==TRUE){
		candin=1
	}else{candin=0}
	return(candin)

}

registerDoMC(cores = 2)
# create dummy expression zscore datafrmae by shuffling expression value across genes
cand_genes_permut_noexp=foreach(i = 1:10000, .combine = rbind) %dopar% {
	
	#zscore_dummy=exp_outlier
	#zscore_dummy$zscore=ave(exp_outlier$zscore,exp_outlier$gene, FUN = sample)

	# select only samples with gendata
	#zscore_dummy=zscore_dummy %>% filter(sample_id %in% metadata_udn_cand$sample_id)
	# merge zscore dummy with RV file

	#RV_dummy=RV_outlier_exp_filt_withsgl_10kb_noZ %>% left_join(zscore_dummy, by=c("gene", "indiv_id"="sample_id"), all.x=T)
	candin=sapply(names(CAND),is_can_in_list_noexp)
	sum(candin)/length(names(CAND))*100
}




is_can_in_list_noexp_real=function(sample){
	if (any(CAND[[sample]]%in%results_exp_out[[sample]][[5]][,2])==TRUE){
		candin=1
	}else{candin=0}
	return(candin)

}
candin_noexp=sapply(names(CAND),is_can_in_list_noexp_real)
real_noexp=sum(candin_noexp)/length(names(CAND))*100

candidates_percent=rbind(melt(data.frame(shuffle_candidate=cand_genes_permut)),melt(data.frame(real_data=real)),melt(data.frame(shuffle_gene=cand_genes_permut_noexp)))
candidates_percent$variable=factor(candidates_percent$variable, levels=c("shuffle_candidate","shuffle_gene","real_data"))

#candidates_percent$filter=c(rep("With expression", 10001), rep("Without expression", 10001))

#candidates_percent_plot=ggplot(candidates_percent, aes(x=filter, y=value, color=variable))+ geom_boxplot() + stat_compare_means()

diagnostic_rate=data.frame(solved=c(7.5),strong_candidate=c(16.7), unsolved=c(86.2))
diagnostic_rate.m=melt(diagnostic_rate)
diagnostic_rate_plot=ggplot(diagnostic_rate.m, aes(x=variable, y=value))+geom_bar(stat="identity", fill="grey", color="black")+
	labs(x="", y="Rate (%)") + scale_x_discrete(labels=c("solved" = "Solved", "strong_candidate" = "Strong Candidate", "unsolved" = "Unsolved"))

candidates_percent_plot=ggplot(candidates_percent, aes(x=variable, y=value))+ geom_boxplot(fill="grey", color="black") +
	labs(x="",y="Percentage of cases for which\ncandidate gene is in final list")+
	scale_x_discrete(labels=c("shuffle_candidate" = "Shuffling Candidate", "shuffle_gene" = "Shuffling genes", "real_data" = "Real data"))

FigureS14.plot=ggarrange(diagnostic_rate_plot, candidates_percent_plot,  
          labels = c("A", "B"),
          nrow = 2)

ggsave('FigureS14.pdf', FigureS14.plot, path=paste("./figures_revision/", sep=""), width=6, height=8)
features=c("stop_gained", "stop_lost", "splice_donor_variant","splice_region_variant","protein_altering_variant","frameshift_variant","inframe_deletion","inframe_insertion", "missense_variant","synonymous_variant","intron_variant","upstream_gene_variant","downstream_gene_variant", "5_prime_UTR_variant", "3_prime_UTR_variant","non_coding_transcript_exon_variant","undefined")


variant_distribution=function(sample){
	print(sample)
	annot_short=rv_outlier_annot %>% filter(indiv_id==sample,gene %in% results_exp_out[[sample]][13][[1]][,1]) %>% filter(as.numeric(as.character(phred_cadd))>=10)

	#annot_short=rv_annotation[rv_annotation$indiv_id==sample,gene %in% results_exp_out[[sample]][13][[1]][,1])]
	#samp_exp=RV_outlier_exp_filt_withsgl_10kb_upstream_cadd_rv %>% filter(sample_id==sample, gene %in% results_exp_out[[sample]][13][[1]][,1])
	#if (nrow(annot_short)!=0){
		annot_short$ensembl_tok=factor(annot_short$ensembl_tok, levels=features)
		annot_short=annot_short %>% select(ensembl_tok) %>% table
		#annot_short=annot_short%>% group_by(ensembl_tok) %>% summarize(count=n()) %>% as.data.frame
		#annot_short=replace(annot_short,is.na(annot_short),"N/A")
		#rownames(annot_short)=annot_short[,1]
		#rown=annot_short[,1]
		#annot_short=as.data.frame(t(annot_short))[-1,]

		#if(length(annot_short)==1){
		#	 annot_short=data.frame(toto=annot_short)
		#	 colnames(annot_short)=rown
		#}
		#samp_exp=samp_exp[,-1]
	#}
	#else{
	#	annot_short=data.frame(matrix(ncol=0,nrow=0))}
	return(annot_short)
}

variant_dist=sapply(cases_withGen,variant_distribution)


variant_dist.m=melt(variant_dist)

getPalette = colorRampPalette(brewer.pal(9, "Set3"))
color.ramp <- color.function(17)

variant_annotation_plot=ggplot(variant_dist.m, aes(x=Var2, y=as.numeric(value), fill=Var1))+ geom_bar(stat="identity") + 
	scale_fill_manual(values= getPalette(17),name="Annotation") + 
	theme(axis.text.x=element_blank())+ 
	labs(x="Samples", y="Number of rare variants")

ggsave('variant_annotation_plot.pdf', variant_annotation_plot, path='figures_revision/', width=12, height=6)


save.image(file=paste0(dir,"/analysis/manuscript/figures_revision/Figure2.in.RData"))




genes_samp_1=RV_withsgl_10kb_annot %>% filter(indiv_id=="RD182", variant_gene_pos2=="upstream", as.numeric(as.character(phred_cadd))>=10) %>% mutate(ensgene=str_extract(gene, "ENSG[0-9]+")) %>% left_join(grch37,by="ensgene") %>% left_join(pheno_annot, by=c("entrez", "symbol")) %>%mutate(hpo_in_input=ifelse(HPO_ID%in% HPO_extended[["RD182"]], 1, 0))  %>% group_by(ensgene) %>% summarize(sum_=sum(hpo_in_input))%>% filter(sum_>0)%>% select(ensgene) %>%pull
genes_samp_2=RV_withsgl_10kb_annot %>% filter(indiv_id=="RD194", variant_gene_pos2=="upstream", as.numeric(as.character(phred_cadd))>=10) %>% mutate(ensgene=str_extract(gene, "ENSG[0-9]+")) %>% left_join(grch37,by="ensgene") %>% left_join(pheno_annot, by=c("entrez", "symbol")) %>%mutate(hpo_in_input=ifelse(HPO_ID%in% HPO_extended[["RD194"]], 1, 0))  %>% group_by(ensgene) %>% summarize(sum_=sum(hpo_in_input))%>% filter(sum_>0)%>% select(ensgene) %>%pull

RV_withsgl_10kb_annot %>% filter(indiv_id=="RD182", variant_gene_pos2=="upstream", as.numeric(as.character(phred_cadd))>=10)%>%  mutate(ensgene=str_extract(gene, "ENSG[0-9]+")) %>% filter(ensgene %in% genes_samp_1) %>% select(ensgene, chr,pos,nvar, ref,alt,phred_cadd,gnomAD_AF_filt_withsgl) %>% group_by(ensgene)%>% summarize(n_=n())%>% filter(n_>1)
RV_withsgl_10kb_annot %>% filter(indiv_id=="RD194", variant_gene_pos2=="upstream", as.numeric(as.character(phred_cadd))>=10)%>%  mutate(ensgene=str_extract(gene, "ENSG[0-9]+")) %>% filter(ensgene %in% genes_samp_2) %>% select(ensgene, chr,pos,nvar, ref,alt,phred_cadd,gnomAD_AF_filt_withsgl) %>% group_by(ensgene)%>% summarize(n_=n())%>% filter(n_>1)

genes_samp_3=RV_withsgl_10kb_annot %>% filter(indiv_id=="RD059", variant_gene_pos2=="upstream", as.numeric(as.character(phred_cadd))>=10) %>% mutate(ensgene=str_extract(gene, "ENSG[0-9]+")) %>% left_join(grch37,by="ensgene") %>% left_join(pheno_annot, by=c("entrez", "symbol")) %>%mutate(hpo_in_input=ifelse(HPO_ID%in% HPO_extended[["RD059"]], 1, 0))  %>% group_by(ensgene) %>% summarize(sum_=sum(hpo_in_input))%>% filter(sum_>0)%>% select(ensgene) %>%pull
genes_samp_4=RV_withsgl_10kb_annot %>% filter(indiv_id=="RD062", variant_gene_pos2=="upstream", as.numeric(as.character(phred_cadd))>=10) %>% mutate(ensgene=str_extract(gene, "ENSG[0-9]+")) %>% left_join(grch37,by="ensgene") %>% left_join(pheno_annot, by=c("entrez", "symbol")) %>%mutate(hpo_in_input=ifelse(HPO_ID%in% HPO_extended[["RD062"]], 1, 0))  %>% group_by(ensgene) %>% summarize(sum_=sum(hpo_in_input))%>% filter(sum_>0)%>% select(ensgene) %>%pull



RV_withsgl_10kb_annot %>% filter(indiv_id=="RD059", variant_gene_pos2=="upstream", as.numeric(as.character(phred_cadd))>=10)%>%  mutate(ensgene=str_extract(gene, "ENSG[0-9]+")) %>% filter(ensgene %in% genes_samp_3) %>% select(ensgene, chr,pos,nvar, ref,alt,phred_cadd,gnomAD_AF_filt_withsgl) %>% group_by(ensgene)%>% summarize(n_=n())%>% filter(n_>1)
RV_withsgl_10kb_annot %>% filter(indiv_id=="RD062", variant_gene_pos2=="upstream", as.numeric(as.character(phred_cadd))>=10)%>%  mutate(ensgene=str_extract(gene, "ENSG[0-9]+")) %>% filter(ensgene %in% genes_samp_4) %>% select(ensgene, chr,pos,nvar, ref,alt,phred_cadd,gnomAD_AF_filt_withsgl) %>% group_by(ensgene)%>% summarize(n_=n())%>% filter(n_>1)
