library(data.table)
library(dplyr)
library(ggplot2)
require(doMC)
require(foreach)
library(argparse)
library(viridis)

parser = ArgumentParser()
parser$add_argument('--datadir', help = 'Local output directory')
parser$add_argument('--ameliedir', help = 'Directory with Amelie scripts')
parser$add_argument('--rdsdir', help = 'Rare disease data directory')
parser$add_argument('--topN', help = 'Number of ASE outliers to keep')

args = parser$parse_args()

data_dir = args$datadir
amelie_dir = args$ameliedir
rds_dir = args$rdsdir
topN = as.numeric(args$topN)

registerDoMC(cores = 20)

# # Read in output of get_ASE_outliers_withGTEx.R
gtex_rds_ase_zscores = fread(paste0(data_dir, 'RDS_GTEX_ASE_zscores.txt'))
# Read in metadata file
metadata = fread(paste0(rds_dir, 'metadata/2018_12_02_Rare_Disease_Metadata.tsv'), data.table=F)
conids = filter(metadata, affected_status == 'Case', in_freeze == 'yes')$sample_id
conids2 = filter(metadata, affected_status == 'Control', in_freeze == 'yes')$sample_id
allGenes = unique(gtex_rds_ase_zscores$GeneName)

# # Read in HPO data - NEW
gene2annot = fread(paste0(data_dir, '18_10_23_ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt'),header=F,sep='\t')
colnames(gene2annot) = c('entrez','Gene_symbol','Pheno','HPO_terms_ID')
parent_terms = fread(paste0(data_dir, 'parent_terms_hpo_2018_10_23.txt'),header=F)
child_terms = fread(paste0(data_dir, 'child_terms_hpo_2018_10_23.txt'),header=F)

get_hpo_genes = function(csample,terms=F) {
  tryCatch({
    sample_data = filter(metadata, sample_id == csample)
    hpo_terms = strsplit(sample_data$HPO_terms_ID, ',')[[1]]
    parent_hpo = unlist(sapply(hpo_terms, function(x) strsplit(filter(parent_terms,V1 ==x)$V2,',')[[1]]))
    child_hpo = unlist(sapply(hpo_terms, function(x) strsplit(filter(child_terms,V1 ==x)$V2,',')[[1]]))
    hpo_terms = c(hpo_terms,parent_hpo,child_hpo)
    names(hpo_terms) = NULL
    if (terms) {
      if (length(hpo_terms) == 0) { return(c(csample,NA)) }
      return(cbind(csample,hpo_terms))
    }
    sample_genes = filter(gene2annot, HPO_terms_ID %in% all_terms)$Gene_symbol
    ensg_genes = unique(filter(gen.mapper, HGNC %in% sample_genes)$ENSG)
    return(cbind(csample, ensg_genes))
  }, error = function(e) {
    print(csample)
    return(c(csample,NA))
  })
}

getCounts <- function(case_outliers,allGenes,varFlag) {
  matchCount = data.frame(ID = character(), NumMatch = numeric(), Category = character())
  for (cid in unique(case_outliers$SampID)) {
    chpo = filter(all_hpo_terms, csample == cid)$hpo_terms
    hgenes = filter(gene2annot,HPO_terms_ID %in% chpo)$Gene_symbol
    ensg_genes = filter(gene_mapper,gene %in% hgenes)$ENSG
    ensg_genes = sapply(ensg_genes, function(x) strsplit(x,'[.]')[[1]][1])
    outlier_hpo = filter(case_outliers,SampID==cid,GeneName %in% ensg_genes)
    matchCount = rbind(matchCount, data.frame(ID = cid, NumMatch = nrow(outlier_hpo), Category = 'Matched'))
  }

  randomCount = foreach(i = 1:100, .combine = rbind) %dopar% {
    rCount = data.frame(ID = character(), NumMatch = numeric(), Category = character())
    for (cid in unique(case_outliers$SampID)) {
      chpo = filter(all_hpo_terms, csample == cid)$hpo_terms
      hgenes = filter(gene2annot,HPO_terms_ID %in% chpo)$Gene_symbol
      ensg_genes = filter(gene_mapper,gene %in% hgenes)$ENSG
      ensg_genes = sapply(ensg_genes, function(x) strsplit(x,'[.]')[[1]][1])
      if (varFlag == 'RV') {
        allGenes = unique(filter(udn_variants,SampID==cid)$GeneName)
        randomGenes = sample(allGenes,topN)
      } else if (varFlag == 'CRV') {
        allGenes = unique(filter(udn_variants,SampID==cid,phred_cadd>=10)$GeneName)
        randomGenes = sample(allGenes,topN)
      } else {
        randomGenes = sample(allGenes,topN)
      }
      outlier_hpo = filter(case_outliers,SampID==cid,GeneName %in% randomGenes)
      rCount = rbind(rCount,data.frame(ID = cid, NumMatch = nrow(outlier_hpo), Category = 'Random'))
    }
    return(rCount)
  }
  return(rbind(matchCount,randomCount))
}

all_hpo_terms = do.call(rbind, lapply(filter(metadata, affected_status == 'Case', in_freeze == 'yes')$sample_id, function(x) get_hpo_genes(x,TRUE)))
all_hpo_terms = as.data.frame(all_hpo_terms)

rds_gtex_single_outliers = gtex_rds_ase_zscores %>% group_by(GeneName) %>%
  mutate(UpRank = rank(-RatioScaled)) %>% mutate(DownRank = rank(RatioScaled)) %>%
  mutate(BothRank = rank(-abs(RatioScaled))) %>% ungroup()

exac_metrics = fread(paste0(data_dir,'exac_gene_metrics.txt'))
gene_mapper = fread(paste0(data_dir,'gencode.v19.mapper.txt'))

colnames(gene_mapper)[7] = 'gene'
exac_metrics = merge(exac_metrics,gene_mapper,by='gene')
exac_metrics$GeneName = sapply(exac_metrics$ENSG, function(x)
  strsplit(x,'[.]')[[1]][1])
exac_metrics = exac_metrics %>% select(syn_z,mis_z,lof_z,pLI,GeneName)
rds_gtex_single_outliers = merge(rds_gtex_single_outliers,exac_metrics,by='GeneName')

# # Just ASE outlier status
print('ASE outlier status')
case_outliers = rds_gtex_single_outliers %>%
  mutate(Source = ifelse(grepl('RD', SampID), 'RDS', 'GTEX')) %>%
  mutate(Status = ifelse(SampID %in% conids, 1, 0)) %>%
  filter(Status == 1) %>% group_by(SampID,GeneName) %>%
  top_n(1,abs(RatioScaled)) %>% ungroup() %>% group_by(SampID) %>%
  top_n(topN,abs(RatioScaled)) %>% ungroup()
allGenes = unique(rds_gtex_single_outliers$GeneName)

hpoMatchesOutliers = getCounts(case_outliers,allGenes,varFlag='None')

# # pLI
print('ASE outlier + pLI')
case_outliers = rds_gtex_single_outliers %>%
  mutate(Source = ifelse(grepl('RD', SampID), 'RDS', 'GTEX')) %>%
  mutate(Status = ifelse(SampID %in% conids, 1, 0)) %>%
  filter(Status == 1) %>% group_by(SampID,GeneName) %>%
  top_n(1,abs(RatioScaled)) %>% ungroup() %>% group_by(SampID) %>%
  filter(pLI>0.9) %>% top_n(topN,abs(RatioScaled)) %>% ungroup()
allGenes = unique(filter(rds_gtex_single_outliers,pLI > 0.9)$GeneName)

hpoMatchesOutliersPli = getCounts(case_outliers,allGenes,varFlag='None')

# # RV
print('ASE outlier + RV')
udn_variants = fread(paste0(data_dir, 'RV_withCADD.txt'))
udn_variants$gene_id = sapply(udn_variants$gene_id,function(x) strsplit(x,'[.]')[[1]][1])
udn_variants = udn_variants %>% select(gene_id,chr,pos,indiv_id,phred_cadd)

case_outliers = rds_gtex_single_outliers %>%
  mutate(Source = ifelse(grepl('RD', SampID), 'RDS', 'GTEX')) %>%
  mutate(Status = ifelse(SampID %in% conids, 1, 0)) %>%
  filter(Status == 1)
colnames(udn_variants)[1] = 'GeneName'
colnames(udn_variants)[4] = 'SampID'

case_outliers = merge(case_outliers,udn_variants,by=c('GeneName','SampID'))
case_rv_outliers = case_outliers %>% group_by(SampID,GeneName) %>%
  top_n(1,abs(RatioScaled)) %>% ungroup() %>% group_by(SampID) %>%
  top_n(topN,abs(RatioScaled)) %>% ungroup()
allGenes = unique(udn_variants$GeneName)

hpoMatchesOutliersRV = getCounts(case_rv_outliers,allGenes,varFlag='RV')

# # conserved RV
print('ASE outlier + conserved RV')
colnames(udn_variants)[1] = 'GeneName'
colnames(udn_variants)[4] = 'SampID'
case_crv_outliers = case_outliers %>% group_by(SampID,GeneName) %>%
  top_n(1,abs(RatioScaled)) %>% ungroup() %>% group_by(SampID) %>%
  filter(phred_cadd >= 10) %>% top_n(topN,abs(RatioScaled)) %>% ungroup()
allGenes = unique(filter(udn_variants,phred_cadd >= 10)$GeneName)

hpoMatchesOutliersCRV = getCounts(case_crv_outliers,allGenes,varFlag='CRV')

rm(udn_variants)
save.image(paste0(data_dir, 'ASE_top', topN, '_HPO_matches.RData'))

