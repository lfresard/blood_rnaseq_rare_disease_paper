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

get_genes_filters_splicing=function(sample){
	
	## Number of rare variants in the sample
	print(sample)
	MAF=0.001
	## splicing outlier genes
	SPLI_OUTLIER=zscores_blood[zscores_blood$variable==sample & zscores_blood$absZ >=2,]
	ind = apply(SPLI_OUTLIER, 1, function(x) all(is.na(x)))
	SPLI_OUTLIER =SPLI_OUTLIER[ !ind, ]
	SPLI_OUTLIER_nb=data.frame(ensgene=unique(SPLI_OUTLIER$gene))
	SPLI_OUTLIER_nb=distinct(merge(SPLI_OUTLIER_nb, grch37, by="ensgene", all.x=T)[c("ensgene", "symbol")])
	SPLI_OUTLIER_nb=distinct(merge(SPLI_OUTLIER_nb, SPLI_OUTLIER,by.x="ensgene",by.y="gene", all.x=T)[c("ensgene", "symbol", "value")])
	SPLI_OUTLIER_nb=SPLI_OUTLIER_nb[complete.cases(SPLI_OUTLIER_nb),]
	if(nrow(SPLI_OUTLIER)!=0){
		SPLI_OUTLIER_nb=aggregate(SPLI_OUTLIER_nb$value, by=list(SPLI_OUTLIER_nb$ensgene,SPLI_OUTLIER_nb$symbol), FUN=max, na.rm=T)
		colnames(SPLI_OUTLIER_nb)=c("ensgene", "symbol", "maxZ")
		SPLI_OUTLIER_nb=SPLI_OUTLIER_nb[complete.cases(SPLI_OUTLIER_nb),]
		SPLI_OUTLIER_nb=SPLI_OUTLIER_nb[order(-abs(SPLI_OUTLIER_nb$maxZ)),] 
	}else{
		SPLI_OUTLIER_nb=NA
	}
	## splicing outlier genes with pLI >=0.9
	if(is.na(SPLI_OUTLIER_nb)==FALSE){
		SPLI_OUTLIER_pLI=merge(SPLI_OUTLIER,exac, by.x="gene", by.y="ensgene", all.x=T)
		SPLI_OUTLIER_pLI=SPLI_OUTLIER_pLI[SPLI_OUTLIER_pLI$pLI >=0.9,]
		SPLI_OUTLIER_pLI_nb=data.frame(ensgene=unique(SPLI_OUTLIER_pLI$gene))
		SPLI_OUTLIER_pLI_nb=distinct(merge(SPLI_OUTLIER_pLI_nb, grch37, by="ensgene", all.x=T)[c("ensgene", "symbol")])
		SPLI_OUTLIER_pLI_nb=distinct(merge(SPLI_OUTLIER_pLI_nb, SPLI_OUTLIER,by.x="ensgene",by.y="gene", all.x=T)[c("ensgene", "symbol", "value")])
		SPLI_OUTLIER_pLI_nb=SPLI_OUTLIER_pLI_nb[complete.cases(SPLI_OUTLIER_pLI_nb),]
		
		SPLI_OUTLIER_pLI_nb=aggregate(SPLI_OUTLIER_pLI_nb$value, by=list(SPLI_OUTLIER_pLI_nb$ensgene,SPLI_OUTLIER_pLI_nb$symbol), FUN=max, na.rm=T)
		colnames(SPLI_OUTLIER_pLI_nb)=c("ensgene", "symbol", "maxZ")
		SPLI_OUTLIER_pLI_nb=SPLI_OUTLIER_pLI_nb[complete.cases(SPLI_OUTLIER_pLI_nb),]
		SPLI_OUTLIER_pLI_nb=SPLI_OUTLIER_pLI_nb[order(-abs(SPLI_OUTLIER_pLI_nb$maxZ)),] 
	}else{SPLI_OUTLIER_pLI_nb=NA}
	## splicing outlier HPO matched gene
	if (anyNA(HPO_extended[[sample]])||is.na(SPLI_OUTLIER_nb)==TRUE){
		SPLI_OUT_HPO_nb=NA 
		}
	else{

		SPLI_OUT_HPO=candidates_zscore_hpo_ext[candidates_zscore_hpo_ext$sample==sample,]
		SPLI_OUT_HPO_nb=data.frame(symbol=unique(SPLI_OUT_HPO$symbol))
		SPLI_OUT_HPO_nb=distinct(merge(SPLI_OUT_HPO_nb, grch37, by="symbol", all.x=T)[c("ensgene", "symbol")])
		SPLI_OUT_HPO_nb=distinct(merge(SPLI_OUT_HPO_nb, SPLI_OUTLIER,by.x="ensgene",by.y="gene", all.x=T)[c("ensgene", "symbol", "value")])
		SPLI_OUT_HPO_nb=SPLI_OUT_HPO_nb[complete.cases(SPLI_OUT_HPO_nb),]
		if(nrow(SPLI_OUT_HPO)!=0){
			SPLI_OUT_HPO_nb=aggregate(abs(SPLI_OUT_HPO_nb$value), by=list(SPLI_OUT_HPO_nb$ensgene,SPLI_OUT_HPO_nb$symbol), FUN=max, na.rm=T)
			colnames(SPLI_OUT_HPO_nb)=c("ensgene", "symbol", "maxZ")
			SPLI_OUT_HPO_nb=SPLI_OUT_HPO_nb[complete.cases(SPLI_OUT_HPO_nb),]
			SPLI_OUT_HPO_nb=SPLI_OUT_HPO_nb[order(-abs(SPLI_OUT_HPO_nb$maxZ)),] 
		}else{SPLI_OUT_HPO_nb=NA}
	}
	## splicing outlier with RV within 20bp
	if (sample %in% cases_withGen && is.na(SPLI_OUTLIER_nb)==FALSE){
		SPLI_OUT_RV=RV_outlier_filt_withsgl[RV_outlier_filt_withsgl$variable==sample & RV_outlier_filt_withsgl$absZ >=2,]
		ind = apply(SPLI_OUT_RV, 1, function(x) all(is.na(x)))
		SPLI_OUT_RV =SPLI_OUT_RV[ !ind, ]
		SPLI_OUT_RV_nb=data.frame(ensgene=unique(SPLI_OUT_RV$gene))
		SPLI_OUT_RV_nb=distinct(merge(SPLI_OUT_RV_nb, grch37, by="ensgene", all.x=T)[c("ensgene", "symbol")])
		SPLI_OUT_RV_nb=distinct(merge(SPLI_OUT_RV_nb, SPLI_OUT_RV,by.x="ensgene",by.y="gene", all.x=T)[c("ensgene", "symbol", "value")])
		SPLI_OUT_RV_nb=SPLI_OUT_RV_nb[complete.cases(SPLI_OUT_RV_nb),]
		if(nrow(SPLI_OUT_RV)!=0){
			SPLI_OUT_RV_nb=aggregate(abs(SPLI_OUT_RV_nb$value), by=list(SPLI_OUT_RV_nb$ensgene,SPLI_OUT_RV_nb$symbol), FUN=max, na.rm=T)
			colnames(SPLI_OUT_RV_nb)=c("ensgene", "symbol", "maxZ")
			SPLI_OUT_RV_nb=SPLI_OUT_RV_nb[complete.cases(SPLI_OUT_RV_nb),]
			SPLI_OUT_RV_nb=SPLI_OUT_RV_nb[order(-abs(SPLI_OUT_RV_nb$maxZ)),] 
		}else{SPLI_OUT_RV_nb=NA}
		
		if(is.na(SPLI_OUT_RV_nb)==F){
			## splicing outlier with RV with CADD>=10 within 20bp
			SPLI_OUT_RV_CADD=SPLI_OUT_RV$gene[as.numeric(as.character(SPLI_OUT_RV$phred_cadd))>=10]
			SPLI_OUT_RV_CADD_nb=data.frame(ensgene=SPLI_OUT_RV_CADD[complete.cases(SPLI_OUT_RV_CADD)])
			SPLI_OUT_RV_CADD_nb=distinct(merge(SPLI_OUT_RV_CADD_nb, grch37, by="ensgene", all.x=T)[c("ensgene", "symbol")])
			SPLI_OUT_RV_CADD_nb=distinct(merge(SPLI_OUT_RV_CADD_nb, SPLI_OUT_RV,by.x="ensgene",by.y="gene", all.x=T)[c("ensgene", "symbol", "value")])
			SPLI_OUT_RV_CADD_nb=SPLI_OUT_RV_CADD_nb[complete.cases(SPLI_OUT_RV_CADD_nb),]
			if(nrow(SPLI_OUT_RV_CADD_nb)!=0){
				SPLI_OUT_RV_CADD_nb=aggregate(abs(SPLI_OUT_RV_CADD_nb$value), by=list(SPLI_OUT_RV_CADD_nb$ensgene,SPLI_OUT_RV_CADD_nb$symbol), FUN=max, na.rm=T)
				colnames(SPLI_OUT_RV_CADD_nb)=c("ensgene", "symbol", "maxZ")
				SPLI_OUT_RV_CADD_nb=SPLI_OUT_RV_CADD_nb[complete.cases(SPLI_OUT_RV_CADD_nb),]
				SPLI_OUT_RV_CADD_nb=SPLI_OUT_RV_CADD_nb[order(-abs(SPLI_OUT_RV_CADD_nb$maxZ)),] 
			}else{SPLI_OUT_RV_CADD_nb=NA}	
		}else{SPLI_OUT_RV_CADD_nb=NA}
	## splicing outlier matching HPO with RV with CADD>=10 within 20bp
	if (anyNA(HPO_extended[[sample]]) || is.na(SPLI_OUT_RV_nb)==T) {
		SPLI_OUT_RV_CADD_HPO_nb=NA 
		}
	else{

		SPLI_OUT_RV_CADD_HPO=candidates_zscore_hpo_RV[candidates_zscore_hpo_RV$sample==sample,]
		SPLI_OUT_RV_CADD_HPO_nb=data.frame(symbol=unique(SPLI_OUT_RV_CADD_HPO$symbol))
		SPLI_OUT_RV_CADD_HPO_nb=distinct(merge(SPLI_OUT_RV_CADD_HPO_nb, grch37, by="symbol", all.x=T)[c("ensgene", "symbol")])
		SPLI_OUT_RV_CADD_HPO_nb=distinct(merge(SPLI_OUT_RV_CADD_HPO_nb, SPLI_OUT_RV,by.x="ensgene",by.y="gene", all.x=T)[c("ensgene", "symbol", "value")])
		SPLI_OUT_RV_CADD_HPO_nb=SPLI_OUT_RV_CADD_HPO_nb[complete.cases(SPLI_OUT_RV_CADD_HPO_nb),]
		if(nrow(SPLI_OUT_RV_CADD_HPO_nb)!=0){
			SPLI_OUT_RV_CADD_HPO_nb=aggregate(abs(SPLI_OUT_RV_CADD_HPO_nb$value), by=list(SPLI_OUT_RV_CADD_HPO_nb$ensgene,SPLI_OUT_RV_CADD_HPO_nb$symbol), FUN=max, na.rm=T)
			colnames(SPLI_OUT_RV_CADD_HPO_nb)=c("ensgene", "symbol", "maxZ")
			SPLI_OUT_RV_CADD_HPO_nb=SPLI_OUT_RV_CADD_HPO_nb[complete.cases(SPLI_OUT_RV_CADD_HPO_nb),]
			SPLI_OUT_RV_CADD_HPO_nb=SPLI_OUT_RV_CADD_HPO_nb[order(-abs(SPLI_OUT_RV_CADD_HPO_nb$maxZ)),] 
		}else{SPLI_OUT_RV_CADD_HPO_nb=NA}
	}
	
	}else{
	SPLI_OUT_RV_nb=NA
	SPLI_OUT_RV_CADD_nb=NA
	SPLI_OUT_RV_CADD_HPO_nb=NA

	}
		dfrow = list(SAMPLE=sample,SPLI_OUTLIER=SPLI_OUTLIER_nb,SPLI_OUTLIER_PLI=SPLI_OUTLIER_pLI_nb, SPLI_OUTLIER_HPO= SPLI_OUT_HPO_nb, SPLI_OUTLIER_RV=SPLI_OUT_RV_nb,SPLI_OUTLIER_RV_CADD=SPLI_OUT_RV_CADD_nb,SPLI_OUTLIER_RV_CADD_HPO=SPLI_OUT_RV_CADD_HPO_nb )

	return(dfrow)
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
gene_to_pheno="18_10_23_ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt"
gene_to_pheno = read_tsv(gene_to_pheno, comment= "#",col_names = FALSE)
colnames(gene_to_pheno)=c("entrez", "symbol", "DiseaseId", "HPO_ID")

pheno_annot="18_10_23_ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt"
pheno_annot=read_tsv(pheno_annot, comment= "#",col_names = FALSE)
colnames(pheno_annot)=c("HPO_ID", "HPO_name", "entrez", "symbol")

# Read in and process alternative HPO files (generated by parse_hpoterms.py)
parent_terms="parent_terms_hpo_2018_10_23.txt"
child_terms="child_terms_hpo_2018_10_23.txt"
altid_terms="altid_terms_hpo_2018_10_23.txt"

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


#Read in splicing zscores
zscores_blood =read_tsv('2018_11_13_splicing_zscores_sorted.txt')

# Read in outliers filtered by RV
RV_outlier_file='./RV/2018_11_13_splicing_zscores_sorted_RV_window.txt'
RV_outlier=as.data.frame(read.table(RV_outlier_file, sep="\t", header=F))
colnames(RV_outlier)=c(colnames(zscores_blood), "chr_var", "start_var", "end_var", "ref_var", "alt_var", "gnomAD_AF", "raw_cadd", "phred_cadd", "chr_gene", "start_gene", "end_gene", "gene_name")

	# make a column with singleton included 
RV_outlier =RV_outlier %>% mutate(gnomAD_AF_withsgl = replace(gnomAD_AF, gnomAD_AF==".", 0))

# format variant data so that only the max observed allele frequency is kept for variant were several frequencies are available
	RV_outlier$gnomAD_AF_filt_withsgl=sapply(RV_outlier$gnomAD_AF_withsgl,process_variant)

#filter for RV with MAF<=0.1%
maf=0.001
RV_outlier_filt_withsgl=RV_outlier %>% filter(!is.na(gnomAD_AF_filt_withsgl)) %>% filter(gnomAD_AF_filt_withsgl <=maf)

	# change gene name format
RV_outlier_filt_withsgl =RV_outlier_filt_withsgl %>% mutate(gene_name=str_extract(gene_name, "ENSG[0-9]+"))


# generate results with RV
potential_can=lapply(cases_withGen,get_hpo_comp_window,RV_outlier_filt_withsgl_annot,metadata, gene_to_pheno,pheno_annot, 10,HPO_extended)
names(potential_can)=cases_withGen

candidates_zscore_hpo_RV=do.call("rbind", potential_can)
candidates_zscore_hpo_RV=candidates_zscore_hpo_RV %>% arrange(sample)%>% left_join(metadata, by=c("sample"="sample_id")) %>% select( symbol,overlap_hpo,    max_z ,sample, indv_id, institution)


# Get list of genes at each filter
res_splicing_genes=mclapply(cases_withGen,get_genes_filters_splicing, mc.cores=10)
names(res_splicing_genes)=cases_withGen
indiv_id=data.frame(STAN_ID=cases_withGen)
indiv_id=indiv_id %>% left_join(metadata, by=c("STAN_ID"="sample_id")) %>% select(STAN_ID,institution_id)



# write down all splicing results in serated files for each sample
for (i in c(1:length(to_print))){
	print(i)
	for (j in c(2:length(names(res_splicing_genes[[to_print[i]]])))) {
		write.table(res_splicing_genes[[to_print[i]]][[j]], file=paste0(dir,indiv_id[i,2],"_",names(res_splicing_genes[[to_print[i]]][j]),"_splicing_filtered_results.tsv"), sep="\t", quote=F, row.names=F)
	}
}

