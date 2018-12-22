
#!/bin/R


# Master directory
dir = Sys.getenv('RARE_DIS_DIR')

(window method)
RV_10kb=paste(dir,"/analysis/outlier_analysis/expression_level/rare_var_protein_coding_genes.txt", sep="")
RV_10kb=as.data.frame(read.table(RV_10kb, sep="\t", header=F))
colnames(RV_10kb)=c("sample_id", "gene",  "variantPos", "gnomadMAF", "raw_cadd", "phred_cadd")

# make a column with singleton included or expluded form the analysis (window method)
RV_10kb =RV_10kb %>% mutate(gnomAD_AF_withsgl = replace(gnomadMAF, gnomadMAF==".", 0))
RV_10kb =RV_10kb %>% mutate(gnomAD_AF_nosgl = replace(gnomadMAF, gnomadMAF==".", NA))


# format variant data so that only the max observed allele frequency is kept for variant were several frequencies are available (window method)
RV_10kb$gnomAD_AF_filt_withsgl=sapply(RV_10kb$gnomAD_AF_withsgl,process_variant)

#filter for RV with MAF<=0.1% (window method)
maf=0.001
	# Including singletons
RV_withsgl_10kb=RV_10kb %>% filter(!is.na(gnomAD_AF_filt_withsgl)) %>% filter(gnomAD_AF_filt_withsgl <=maf)


# get number of genes passing different criterion for sample of interest.
rars2_RD059.df=data.frame(
	RV_10kb=RV_withsgl_10kb%>% filter(sample_id=="RD059")%>% mutate(gene=str_extract(gene, "ENSG[0-9]+"))%>% left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=0, pLI>=0)  %>% select(gene) %>% mutate(gene=str_extract(gene, "ENSG[0-9]+"))%>% filter (! duplicated(gene)) %>% nrow,
	RV_10kb_CADD=RV_withsgl_10kb%>% filter(sample_id=="RD059")%>% mutate(gene=str_extract(gene, "ENSG[0-9]+"))%>% left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=10, pLI>=0)  %>% select(gene) %>% mutate(gene=str_extract(gene, "ENSG[0-9]+"))%>% filter (! duplicated(gene)) %>% nrow,
	RV_10kb_HPO=RV_withsgl_10kb%>% filter(sample_id=="RD059")%>% mutate(gene=str_extract(gene, "ENSG[0-9]+")) %>%  left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct %>% filter(pLI>=0) %>%	left_join(grch37, by=c("gene"="ensgene")) %>%	select(gene, symbol, entrez) %>% 	left_join(gene_to_pheno, by=c("entrez","symbol"))%>% distinct %>% 	left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>%		mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPO_extended[["RD059"]], 1, 0)) %>%select(gene, symbol,  hpo_in_input) %>% distinct %>%group_by(gene,symbol)%>% summarise(overlap_hpo=sum(hpo_in_input)) %>%  		filter(overlap_hpo>=1) %>% ungroup %>% unique%>% filter (! duplicated(gene)) %>% nrow,
	RV_10kb_CADD_HPO=RV_withsgl_10kb%>% filter(sample_id=="RD059")%>% mutate(gene=str_extract(gene, "ENSG[0-9]+")) %>% filter(as.numeric(as.character(phred_cadd))>=10)%>% left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct %>% filter(pLI>=0) %>%	left_join(grch37, by=c("gene"="ensgene")) %>%	select(gene, symbol, entrez) %>% 	left_join(gene_to_pheno, by=c("entrez","symbol"))%>% distinct %>% 	left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>%		mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPO_extended[["RD059"]], 1, 0)) %>%select(gene, symbol, hpo_in_input) %>% distinct %>%group_by(gene,symbol)%>% summarise(overlap_hpo=sum(hpo_in_input)) %>%  		filter(overlap_hpo>=1) %>% ungroup %>% unique%>% filter (! duplicated(gene)) %>% nrow,
	hpo_match=gene_to_pheno %>% left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>% mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPO_extended[["RD059"]], 1, 0)) %>% filter(hpo_in_input==1) %>% select(entrez) %>% unique %>% nrow,
	exp_outlier=exp_outlier %>%filter(sample_id=="RD059",abs(zscore)>=2) %>% unique %>%nrow,
	exp_out_pLI=exp_outlier %>%filter(sample_id=="RD059",abs(zscore)>=2) %>% unique%>% left_join(exac,by=c("gene"="ensgene"), all.x=T)%>% distinct %>% filter(pLI>=0.9) %>% select(gene, sample_id, zscore) %>% unique  %>%nrow,
	exp_out_river=exp_outlier %>%filter(sample_id=="RD059",abs(zscore)>=2, gene %in% gene_river_list[["RD059"]]) %>% unique %>% filter (! duplicated(gene))%>%nrow,
	exp_out_ase=exp_outlier %>%filter(sample_id=="RD059",abs(zscore)>=2, gene %in% gene_ase_list[["RD059"]]) %>% unique %>% filter (! duplicated(gene))%>%nrow,
	exp_out_hpo=exp_outlier %>% filter(sample_id=="RD059")  %>% filter(abs(zscore)>=2) %>% left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct %>% filter(pLI>=0) %>%	left_join(grch37, by=c("gene"="ensgene")) %>% select(gene, symbol, entrez, zscore) %>% 	left_join(gene_to_pheno, by=c("entrez","symbol"))%>% distinct %>% 	left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>%	mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPO_extended[["RD059"]], 1, 0)) %>%select(gene, symbol, zscore, hpo_in_input) %>% distinct %>%group_by(gene,symbol, zscore)%>% summarise(overlap_hpo=sum(hpo_in_input)) %>% filter(overlap_hpo>=1) %>% ungroup %>% unique%>% filter (! duplicated(gene)) %>% nrow,
	exp_out_RV=RV_outlier_exp_filt_withsgl_10kb %>% filter(sample_id=="RD059")%>% filter(abs(expressionZScore)>=2) %>% left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=0, pLI>=0)  %>% select(gene,expressionZScore) %>%unique%>% filter (! duplicated(gene)) %>% nrow,
	exp_out_RV_CADD=RV_outlier_exp_filt_withsgl_10kb %>% filter(sample_id=="RD059")%>% filter(abs(expressionZScore)>=2) %>% filter(as.numeric(as.character(phred_cadd))>=10)%>% left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=0, pLI>=0)  %>% select(gene,expressionZScore) %>%unique%>% filter (! duplicated(gene)) %>% nrow,
	exp_out_RV_CADD_HPO=RV_outlier_exp_filt_withsgl_10kb%>% filter(sample_id=="RD059") %>% filter(abs(expressionZScore)>=2) %>% filter(as.numeric(as.character(phred_cadd))>=10)%>% left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct %>% filter(pLI>=0) %>%	left_join(grch37, by=c("gene"="ensgene")) %>%	select(gene, symbol, entrez, expressionZScore) %>% 	left_join(gene_to_pheno, by=c("entrez","symbol"))%>% distinct %>% 	left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>%		mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPO_extended[["RD059"]], 1, 0)) %>%select(gene, symbol, expressionZScore, hpo_in_input) %>% distinct %>%group_by(gene,symbol, expressionZScore)%>% summarise(overlap_hpo=sum(hpo_in_input)) %>%  		filter(overlap_hpo>=1) %>% ungroup %>% unique%>% filter (! duplicated(gene)) %>% nrow)
rars2_RD059.df.m=melt( rars2_RD059.df)
rars2_RD059.df.m$category=c("Rare variant", rep("Rare variant and filters",3),"HPO", "Expression outlier", rep("Expression and filters",7))

rars2_RD059.df.m$category=factor(rars2_RD059.df.m$category, levels=c("Expression and filters","Expression outlier","HPO","Rare variant and filters","Rare variant"))


# Look if causal gene is in the list of filtered genes
rars2_RD059.causalgene.df=data.frame(
	RV_10kb=RV_withsgl_10kb%>% filter(sample_id=="RD059")%>% mutate(gene=str_extract(gene, "ENSG[0-9]+"))%>% left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=0, pLI>=0)  %>% select(gene) %>% mutate(gene=str_extract(gene, "ENSG[0-9]+"))%>% filter (! duplicated(gene)) %>%mutate(is_causal=ifelse("ENSG00000146282" %in% gene, 1, 0)) %>% select(is_causal) %>% distinct %>% pull,
	RV_10kb_CADD=RV_withsgl_10kb%>% filter(sample_id=="RD059")%>% mutate(gene=str_extract(gene, "ENSG[0-9]+"))%>% left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=10, pLI>=0)  %>% select(gene) %>% mutate(gene=str_extract(gene, "ENSG[0-9]+"))%>% filter (! duplicated(gene))  %>%mutate(is_causal=ifelse("ENSG00000146282" %in% gene, 1, 0)) %>% select(is_causal) %>% distinct %>% pull,
	RV_10kb_HPO=RV_withsgl_10kb%>% filter(sample_id=="RD059")%>% mutate(gene=str_extract(gene, "ENSG[0-9]+")) %>%  left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct %>% filter(pLI>=0) %>%	left_join(grch37, by=c("gene"="ensgene")) %>%	select(gene, symbol, entrez) %>% 	left_join(gene_to_pheno, by=c("entrez", "symbol"))%>% distinct %>% 	left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>%		mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPO_extended[["RD059"]], 1, 0)) %>%select(gene, symbol,  hpo_in_input) %>% distinct %>%group_by(gene,symbol)%>% summarise(overlap_hpo=sum(hpo_in_input)) %>%  		filter(overlap_hpo>=1) %>% ungroup %>% unique%>% filter (! duplicated(gene)) %>%mutate(is_causal=ifelse("ENSG00000146282" %in% gene, 1, 0)) %>% select(is_causal) %>% distinct %>% pull,
	RV_10kb_CADD_HPO=RV_withsgl_10kb%>% filter(sample_id=="RD059")%>% mutate(gene=str_extract(gene, "ENSG[0-9]+")) %>% filter(as.numeric(as.character(phred_cadd))>=10)%>% left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct %>% filter(pLI>=0) %>%	left_join(grch37, by=c("gene"="ensgene")) %>%	select(gene, symbol, entrez) %>% 	left_join(gene_to_pheno, by=c("entrez","symbol"))%>% distinct %>% 	left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>%		mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPO_extended[["RD059"]], 1, 0)) %>%select(gene, symbol, hpo_in_input) %>% distinct %>%group_by(gene,symbol)%>% summarise(overlap_hpo=sum(hpo_in_input)) %>%  		filter(overlap_hpo>=1) %>% ungroup %>% unique%>% filter (! duplicated(gene)) %>%mutate(is_causal=ifelse("ENSG00000146282" %in% gene, 1, 0)) %>% select(is_causal) %>% distinct %>% pull,
	hpo_match=gene_to_pheno %>% left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>% mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPO_extended[["RD059"]], 1, 0)) %>% filter(hpo_in_input==1) %>% select(symbol) %>% unique %>% mutate(is_causal=ifelse("RARS2" %in% symbol, 1, 0)) %>% select(is_causal) %>% distinct %>% pull,
	exp_outlier=exp_outlier %>%filter(sample_id=="RD059",abs(zscore)>=2) %>% unique %>%mutate(is_causal=ifelse("ENSG00000146282" %in% gene, 1, 0)) %>% select(is_causal) %>% distinct %>% pull,
	exp_out_pLI=exp_outlier %>%filter(sample_id=="RD059",abs(zscore)>=2) %>% unique%>% left_join(exac,by=c("gene"="ensgene"), all.x=T)%>% distinct %>% filter(pLI>=0.9) %>% select(gene, sample_id, zscore) %>% unique %>%mutate(is_causal=ifelse("ENSG00000146282" %in% gene, 1, 0)) %>% select(is_causal) %>% distinct %>% pull,
	exp_out_river=exp_outlier %>%filter(sample_id=="RD059",abs(zscore)>=2,gene %in% gene_river_list[["RD059"]]) %>% unique %>%mutate(is_causal=ifelse("ENSG00000146282" %in% gene, 1, 0)) %>% select(is_causal) %>% distinct %>% pull,
	exp_out_ase=exp_outlier %>%filter(sample_id=="RD059",abs(zscore)>=2,gene %in% gene_ase_list[["RD059"]]) %>% unique %>%mutate(is_causal=ifelse("ENSG00000146282" %in% gene, 1, 0)) %>% select(is_causal) %>% distinct %>% pull,
	exp_out_hpo=exp_outlier %>% filter(sample_id=="RD059")  %>% filter(abs(zscore)>=1.5) %>% left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct %>% filter(pLI>=0) %>%	left_join(grch37, by=c("gene"="ensgene")) %>% select(gene, symbol, entrez, zscore) %>% 	left_join(gene_to_pheno, by=c("entrez","symbol"))%>% distinct %>% 	left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>%	mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPO_extended[["RD059"]], 1, 0)) %>%select(gene, symbol, zscore, hpo_in_input) %>% distinct %>%group_by(gene,symbol, zscore)%>% summarise(overlap_hpo=sum(hpo_in_input)) %>% filter(overlap_hpo>=1) %>% ungroup %>% unique%>% filter (! duplicated(gene)) %>%mutate(is_causal=ifelse("ENSG00000146282" %in% gene, 1, 0)) %>% select(is_causal) %>% distinct %>% pull,
	exp_out_RV=RV_outlier_exp_filt_withsgl_10kb %>% filter(sample_id=="RD059")%>% filter(abs(expressionZScore)>=2) %>% left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=0, pLI>=0)  %>% select(gene,expressionZScore) %>%unique%>% filter (! duplicated(gene)) %>%mutate(is_causal=ifelse("ENSG00000146282" %in% gene, 1, 0)) %>% select(is_causal) %>% distinct %>% pull,
	exp_out_RV_CADD=RV_outlier_exp_filt_withsgl_10kb %>% filter(sample_id=="RD059")%>% filter(abs(expressionZScore)>=2) %>% filter(as.numeric(as.character(phred_cadd))>=10)%>% left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=0, pLI>=0)  %>% select(gene,expressionZScore) %>%unique%>% filter (! duplicated(gene)) %>%mutate(is_causal=ifelse("ENSG00000146282" %in% gene, 1, 0)) %>% select(is_causal) %>% distinct %>% pull,
	exp_out_RV_CADD_HPO=0)#RV_outlier_exp_filt_withsgl_10kb%>% filter(sample_id=="RD059") %>% filter(abs(expressionZScore)>=2) %>% filter(as.numeric(as.character(phred_cadd))>=10)%>% left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct %>% filter(pLI>=0) %>%	left_join(grch37, by=c("gene"="ensgene")) %>%	select(gene, symbol, entrez, expressionZScore) %>% 	left_join(gene_to_pheno, by=c("entrez","symbol"))%>% distinct %>% 	left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>%		mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPO_extended[["RD059"]], 1, 0)) %>%select(gene, symbol, expressionZScore, hpo_in_input) %>% distinct %>%group_by(gene,symbol, expressionZScore)%>% summarise(overlap_hpo=sum(hpo_in_input)) %>%  		filter(overlap_hpo>=1) %>% ungroup %>% unique%>% filter (! duplicated(gene)) %>%mutate(is_causal=ifelse("ENSG00000146282" %in% gene, 1, 0)) %>% select(is_causal) %>% distinct %>% pull)

rars2_RD059.causalgene.df.m=melt(rars2_RD059.causalgene.df)
rars2_RD059.df.m$is_causal_in=as.factor(rars2_RD059.causalgene.df.m$value)


# Save data
save.image(file = paste(dir,"/data/FigureS15.in.RData",sep=""))
