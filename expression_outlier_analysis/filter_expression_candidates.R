#expression outlier data filters
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(parallel)
library(annotables)
library(cowplot)
library(purrr)
library(broom)
library(reshape2)

### Functions

# format variant data so that only the max observed allele frequency is kept for variant where several frequencies are available 
process_variant=function(variant_freq){
	variant_freq=as.character(variant_freq)
	if (grepl(",", variant_freq)) {
		variant_freq=max(as.numeric(unlist(strsplit(variant_freq,","))))
		return(as.numeric(variant_freq))
	}
	else{ return(as.numeric(variant_freq))}
}



# generates list of HPO alternative IDs from the input files
get_hpo_alt_id=function(myHPO_ID, altid_terms){
	if(altid_terms$ALT_ID[altid_terms$HPO_ID==myHPO_ID]!=""){
		alt_id=altid_terms %>% filter(HPO_ID==myHPO_ID) %>% mutate(ALT_ID=as.character(ALT_ID)) %>%pull
		hpos_alt=unlist(strsplit(alt_id, ","))
		hpos=c(myHPO_ID,hpos_alt)
	}else{hpos=myHPO_ID}
	return(hpos)

}


# generates HPO term list for each sample
get_list_hpo=function(sample, metadata){
	hpos=metadata %>% filter(sample_id==sample) %>% select(HPO_terms_ID)%>% pull
	hpos=unlist(strsplit(hpos,","))
	return(hpos)
}

# function that look for all HPO ters from a given samples and expend the list with parent, child and alt terms
process_sample_ids_hpo=function(sample, HPO_list){
	all_hpos=lapply(HPO_list[[sample]],process_sample_hpo)
	all_hpos=Reduce(union, all_hpos)
	return(all_hpos)
}

get_genes_filters=function(sample){

	## expression outliers
	EXP_SAMPLE=exp_outlier[exp_outlier$sample_id==sample & exp_outlier$zscore <=-2,] # underexpression outliers only
	EXP_SAMPLE$gene=str_extract(EXP_SAMPLE$gene, "ENSG[0-9]+")
	EXP_OUT=data.frame(ensgene=unique(EXP_SAMPLE$gene))
	EXP_OUT=distinct(merge(EXP_OUT, grch37,by="ensgene", all.x=T)[c("ensgene", "symbol")])
	EXP_OUT=distinct(merge(EXP_OUT, EXP_SAMPLE,by.x="ensgene",by.y="gene", all.x=T)[c("ensgene", "symbol", "zscore")])
	EXP_OUT=EXP_OUT[order(-abs(EXP_OUT$zscore)),] 
	
	## Number of expression outliers with pLI>0.9
	exp_out_pLI=EXP_SAMPLE
	exp_out_pLI=distinct(merge(exp_out_pLI,exac,by.x="gene", by.y="ensgene", all.x=T))
	exp_out_pLI=data.frame(ensgene=unique(exp_out_pLI[which(exp_out_pLI$pLI >=0.9),c("gene")]))
	exp_out_pLI=distinct(merge(exp_out_pLI, grch37,by="ensgene", all.x=T)[c("ensgene", "symbol")])
	exp_out_pLI=distinct(merge(exp_out_pLI, EXP_SAMPLE,by.x="ensgene",by.y="gene", all.x=T)[c("ensgene", "symbol", "zscore")])
	exp_out_pLI=exp_out_pLI[order(-abs(exp_out_pLI$zscore)),] 

	## EXP outliers with HPO match
	if (anyNA(HPO_extended[[sample]])){
		exp_out_hpo=NA 
	}
	else{
		exp_out_hpo=distinct(merge(EXP_SAMPLE,grch37, by.x="gene", by.y="ensgene", all.x=T)[,c("gene", "symbol", "entrez", "zscore")])
		exp_out_hpo=merge(exp_out_hpo,pheno_annot, by=c("entrez", "symbol"), all.x=T)
		exp_out_hpo$hpo_in_input=ifelse(exp_out_hpo$HPO_ID %in% HPO_extended[[sample]], 1, 0) # annotate rows for which HPO match
		exp_out_hpo=unique(exp_out_hpo[,c("gene", "symbol",  "hpo_in_input")])
		exp_out_hpo=aggregate(exp_out_hpo$hpo_in_input, by=list(exp_out_hpo$gene), FUN=sum)
		colnames(exp_out_hpo)=c("gene", "hpo_match")
		exp_out_hpo=data.frame(ensgene=unique(exp_out_hpo$gene[exp_out_hpo$hpo_match>0]))
		exp_out_hpo=distinct(merge(exp_out_hpo, grch37,by="ensgene", all.x=T)[c("ensgene", "symbol")])
		exp_out_hpo=distinct(merge(exp_out_hpo, EXP_SAMPLE,by.x="ensgene",by.y="gene", all.x=T)[c("ensgene", "symbol", "zscore")])
		exp_out_hpo=exp_out_hpo[order(-abs(exp_out_hpo$zscore)),] 

	}
	
	## expression outliers with RV 

	if (sample %in% cases_withGen){
		OUT_RV_SAMPLE=RV_outlier_exp_filt_withsgl_10kb[RV_outlier_exp_filt_withsgl_10kb$sample_id==sample & RV_outlier_exp_filt_withsgl_10kb$expressionZScore <=-2,] #underexpression only
		exp_out_RV=data.frame(ensgene=unique(OUT_RV_SAMPLE$gene))
		exp_out_RV=distinct(merge(exp_out_RV, grch37,by="ensgene", all.x=T)[c("ensgene", "symbol")])
		exp_out_RV=distinct(merge(exp_out_RV, OUT_RV_SAMPLE,by.x="ensgene",by.y="gene", all.x=T)[c("ensgene", "symbol", "expressionZScore")])
		exp_out_RV=exp_out_RV[order(-abs(exp_out_RV$expressionZScore)),] 

	
		## expression outliers with RV with CADD>10
		
		exp_out_RV_CADD=data.frame(ensgene=unique(OUT_RV_SAMPLE$gene[which(as.numeric(as.character(OUT_RV_SAMPLE$phred_cadd))>=10)]))
		exp_out_RV_CADD=distinct(merge(exp_out_RV_CADD, grch37,by="ensgene", all.x=T)[c("ensgene", "symbol")])
		exp_out_RV_CADD=distinct(merge(exp_out_RV_CADD, OUT_RV_SAMPLE,by.x="ensgene",by.y="gene", all.x=T)[c("ensgene", "symbol", "expressionZScore")])
		exp_out_RV_CADD=exp_out_RV_CADD[order(-abs(exp_out_RV_CADD$expressionZScore)),] 

	## Num
			if (anyNA(HPO_extended[[sample]]) || nrow(OUT_RV_SAMPLE)==0){
			exp_out_RV_CADD_HPO=NA 
		}	
		else{
			exp_out_RV_CADD_HPO=distinct(merge(OUT_RV_SAMPLE,grch37, by.x="gene", by.y="ensgene", all.x=T)[,c("gene", "symbol", "entrez", "expressionZScore")])
			exp_out_RV_CADD_HPO=merge(exp_out_RV_CADD_HPO,pheno_annot, by=c("entrez", "symbol"), all.x=T)
			exp_out_RV_CADD_HPO$hpo_in_input=ifelse(exp_out_RV_CADD_HPO$HPO_ID %in% HPO_extended[[sample]], 1, 0) # annotate rows for which HPO match
			exp_out_RV_CADD_HPO=unique(exp_out_RV_CADD_HPO[,c("gene", "symbol",  "hpo_in_input")])
			exp_out_RV_CADD_HPO=aggregate(exp_out_RV_CADD_HPO$hpo_in_input, by=list(exp_out_RV_CADD_HPO$gene), FUN=sum)
			colnames(exp_out_RV_CADD_HPO)=c("gene", "hpo_match")
			exp_out_RV_CADD_HPO=data.frame(ensgene=unique(exp_out_RV_CADD_HPO$gene[exp_out_RV_CADD_HPO$hpo_match>0]))
			exp_out_RV_CADD_HPO=distinct(merge(exp_out_RV_CADD_HPO, grch37,by="ensgene", all.x=T)[c("ensgene", "symbol")])
			exp_out_RV_CADD_HPO=distinct(merge(exp_out_RV_CADD_HPO, OUT_RV_SAMPLE,by.x="ensgene",by.y="gene", all.x=T)[c("ensgene", "symbol", "expressionZScore")])
			exp_out_RV_CADD_HPO=exp_out_RV_CADD_HPO[order(-abs(exp_out_RV_CADD_HPO$expressionZScore)),] 
	
		}
		
		}else{
			exp_out_RV=NA
			exp_out_RV_CADD=NA
			exp_out_RV_CADD_HPO=NA
		}


	sample_res=list(SAMPLE=sampleEXP_OUTLIER=EXP_OUT,EXP_OUTLIER_PLI=exp_out_pLI, EXP_OUTLIER_HPO= exp_out_hpo, EXP_OUTLIER_RV=exp_out_RV,EXP_OUTLIER_RV_CADD=exp_out_RV_CADD,EXP_OUTLIER_RV_CADD_HPO=exp_out_RV_CADD_HPO )
	
	return(sample_res)
}

### MAIN

# Get metadata and filter for samples to analyze
metadata = "2018_12_02_Rare_Disease_Metadata.tsv"
metadata = read_tsv(metadata) %>% filter(in_freeze=="yes") %>% filter(source_of_RNA=="Blood",sequencing_status=="PASSED")
metadata = read_tsv(metadata) %>% filter(sequencing_status=="PASSED", is_RD=="yes",source_of_RNA=="Blood")
sample_with_data=metadata %>% filter(variant_data=="exome" | variant_data=="genome")%>% select(sample_id) %>% pull


# read in exac data
exac = as.data.frame(read.table('forweb_cleaned_exac_r03_march16_z_data_pLI.txt', sep = '\t', stringsAsFactors = F, header = T))
	# select gene namem synonymous-zscore, missense-zscore, loss-of-function-zscore and pLI
exac=exac %>% inner_join(grch37, by=c("gene"="symbol")) %>% select(ensgene,syn_z, mis_z, lof_z, pLI) 

# Read in and process HPO data
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

# Generate extended list of HPO terms for each sample
HPO_extended=lapply(cases, process_sample_ids_hpo, HPO_list=HPOs)
names(HPO_extended)=cases

# Read in expression outlier data
exp_outlier="outliers_zscore_pair_spline_07dec18.txt"
exp_outlier=read_tsv(exp_outlier)

# Read in expression outlier with RV
RV_outlier_exp_10kb="outliers_rare_var_combined_10kb_07dec18.txt"
RV_outlier_exp_10kb=read_tsv(RV_outlier_exp_10kb, col_names=FALSE)
	# add column names
colnames(RV_outlier_exp_10kb)=c("sample_id", "gene", "expressionZScore", "variantPos", "gnomadMAF", "raw_cadd", "phred_cadd", "variant_gene_pos")
	# join data to grch37 to get start and end of gene
RV_outlier_exp_10kb=RV_outlier_exp_10kb %>% left_join(grch37, by=c("gene"="ensgene"), all.x=T)%>%select(sample_id,gene,expressionZScore,gnomadMAF,phred_cadd,variant_gene_pos,variantPos, start,end)
	# calculate distance of variant to the gene
RV_outlier_exp_10kb=RV_outlier_exp_10kb%>% mutate(distance=ifelse((variantPos<=start&variantPos<=end)|(variantPos>=start&variantPos>=end),start-variantPos,0))


# make a column with singleton included
RV_outlier_exp_10kb =RV_outlier_exp_10kb %>% mutate(gnomAD_AF_withsgl = replace(gnomadMAF, gnomadMAF==".", 0))


# format variant data so that only the max observed allele frequency is kept for variant where several frequencies are available 
RV_outlier_exp_10kb$gnomAD_AF_filt_withsgl=sapply(RV_outlier_exp_10kb$gnomAD_AF_withsgl,process_variant)


#filter for RV with MAF<=0.1%
maf=0.001
RV_outlier_exp_filt_withsgl_10kb=RV_outlier_exp_10kb %>% filter(!is.na(gnomAD_AF_filt_withsgl)) %>% filter(gnomAD_AF_filt_withsgl <=maf)


# Get list of genes for each filter step
results_exp_out=mclapply(cases_withGen, get_genes_filters, mc.cores=5)
names(results_exp_out)=cases_withGen
indiv_id=data.frame(STAN_ID=cases_withGen)
indiv_id=indiv_id %>% left_join(metadata, by=c("STAN_ID"="sample_id"))  %>% select(STAN_ID,institution_id)

# write down all splicing results in serated files for each sample, each sheet being one filter --> does not work
for (i in c(1:length(cases_withGen))){
	print(i)
	for (j in c(2:length(names(results_exp_out[[i]])))) {
		write.table(results_exp_out[[i]][[j]], file=paste0(indiv_id[i,2],"_",names(results_exp_out[[i]][j]),"_expression_filtered_results.tsv"), sep="\t", quote=F, row.names=F)
	}
}

	if("EXP_OUTLIER_RV_CADD_HPO" %in% names(results_exp_out[[i]])){
	write.table(results_exp_out[[i]][["EXP_OUTLIER_RV_CADD_HPO"]], file=paste0(indiv_id[i,2],"_expression_filtered_results.tsv"), sep="\t", quote=F, row.names=F)
	}else{next}


