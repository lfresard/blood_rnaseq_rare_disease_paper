library(data.table)
library(dplyr)
library(ggplot2)
require(doMC)
require(foreach)

registerDoMC(cores = 20)

data_dir = '/users/nferraro/data/rds_data/'

# Read in output of get_ASE_outliers_withGTEx.R
gtex_rds_ase_zscores = fread(paste0(data_dir, 'RDS_GTEX_ASE_zscores.txt'))
# Read in metadata file
metadata = fread(paste0(data_dir, '2018_06_12_Rare_Disease_Metadata.tsv'), data.table=F)
conids = filter(metadata, affected_status == 'Case', in_freeze == 'yes')$sample_id
conids2 = filter(metadata, affected_status == 'Control', in_freeze == 'yes')$sample_id
allGenes = unique(fread(paste0(data_dir,'rds_gtex_ase_genes.txt'))$GeneName)

# Read in HPO data
pheno_annot = fread(paste0(data_dir, 'phenotype_annotation.tsv'))
gene2pheno = fread(paste0(data_dir, 'gene_to_phenotype.tsv'))
colnames(pheno_annot)[5] = 'DiseaseId'
gene2annot = merge(pheno_annot, gene2pheno, by='DiseaseId')
colnames(gene2annot)[12] = 'Gene_symbol'

get_hpo_genes = function(csample,terms=F) {
  tryCatch({
    sample_data = filter(metadata, sample_id == csample)
    hpo_terms = strsplit(sample_data$HPO_terms_ID, ',')[[1]]
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

all_hpo_terms = do.call(rbind, lapply(filter(metadata, affected_status == 'Case', in_freeze == 'yes')$sample_id, function(x) get_hpo_genes(x,TRUE)))
all_hpo_terms = as.data.frame(all_hpo_terms)

case_outliers = gtex_rds_ase_zscores %>%
  mutate(Source = ifelse(grepl('RD', SampID), 'RDS', 'GTEX')) %>%
  mutate(Status = ifelse(SampID %in% conids, 1, 0)) %>%
  filter(Status == 1, UpRank == 1 | DownRank == 1) %>%
  filter(abs(RatioScaled) > 2)

case_outlier_number = as.data.frame(table(case_outliers$SampID))
colnames(case_outlier_number)[1] = 'sample_id'

unsolved_cases = filter(metadata, affected_status == 'Case', resolved_case == 'no')$sample_id
unsolved_outliers = case_outliers %>% filter(SampID %in% unsolved_cases)

amelie_output = c('ID', 'Genes', 'HPO')
for (cid in unique(unsolved_outliers$SampID)) {
  cgenes = paste(filter(unsolved_outliers, SampID == cid)$GeneName,collapse=',')
  chpo = paste(filter(all_hpo_terms, csample == cid)$hpo_terms,collapse=',')
  amelie_output = rbind(amelie_output, c(cid,cgenes,chpo))
}
colnames(amelie_output) = amelie_output[1,]
amelie_output = amelie_output[-1,]
# Write out output needed for Amelie - used below
write.table(filter(as.data.frame(amelie_output), HPO != ''), file=paste0(data_dir, 'amelie_permutations/RDScases_aseOutliers_amelie.txt'),
            sep='\t', quote=F,row.names=F,col.names=F)

outfile = paste0(data_dir, 'amelie_permutations/RDScases_aseOutliers_amelie.txt')
scommand = paste0('bash /users/nferraro/projects/rds_ase_analysis/scripts/call_amelie.sh ', outfile, ' 0')
system(scommand,ignore.stdout=T)

allRDFiles = dir('/users/nferraro/projects/rds_ase_analysis/scripts/', '*_amelie_rank.txt',full=T)

getAllScores <- function(rdf) {
  tryCatch({
    adata = fread(rdf,header=F, sep='\t')
    cid = strsplit(strsplit(rdf,'/')[[1]][8], '_')[[1]][1] #Index (8) is dependent on file path 
    topScores = unlist(sapply(adata$V2, function(x)
      substr(strsplit(x,',')[[1]][1], start=3, stop=8)))
    adata$TopScore = topScores
    adata$ID = rep(cid,nrow(adata))
    return(adata %>% select(V1,TopScore,ID))
  }, error = function(e) {
    return(c(NA,0,NA))
  })
}

allScoreData = do.call(rbind, lapply(allRDFiles, function(x) getAllScores(x)))
colnames(allScoreData) = c('Gene', 'TopScore', 'ID')
allScoreData$TopScore = as.numeric(allScoreData$TopScore)
NOver50 = nrow(filter(allScoreData, TopScore > 50))
print('Number of genes with scores over 50 in matched cases:')
print(NOver50)

## Permute HPO labels over 100 iterations (calls external amelie script)
allHpoIds = as.character(unique(all_hpo_terms$csample))
countOver50 = 0
numOver50 = c()
allRandomScores = foreach(i = 1:100, .combine = rbind) %dopar% {
  amelie_output = c('ID', 'Genes', 'HPO')
  for (cid in unique(unsolved_outliers$SampID)) {
    ngenes = length(filter(unsolved_outliers, SampID == cid)$GeneName)
    rgenes = sample(allGenes,ngenes)
    cgenes = paste(rgenes,collapse=',')
    chpo = paste(filter(all_hpo_terms, csample == cid)$hpo_terms,collapse=',')
    amelie_output = rbind(amelie_output, c(cid,cgenes,chpo))
  }
  colnames(amelie_output) = amelie_output[1,]
  amelie_output = amelie_output[-1,]
  write.table(filter(as.data.frame(amelie_output), HPO != ''), file=paste0(data_dir, 'amelie_permutations/RDSrandom', i, '_aseOutliers_amelie.txt'),
              sep='\t', quote=F,row.names=F,col.names=F)
  outfile = paste0(data_dir, 'amelie_permutations/RDSrandom', i, '_aseOutliers_amelie.txt')
  scommand = paste0('bash /users/nferraro/projects/rds_ase_analysis/scripts/call_amelie.sh ', outfile, ' ', i)
  system(scommand,ignore.stdout=T,ignore.stderr=T)
  allRDFiles = dir('/users/nferraro/projects/rds_ase_analysis/scripts/', paste0(i,'RD*'),full=T)
  allScoreData = do.call(rbind, lapply(allRDFiles, function(x) getAllScores(x)))
  colnames(allScoreData) = c('Gene', 'TopScore', 'ID')
  allScoreData = as.data.frame(allScoreData)
  allScoreData$TopScore = as.numeric(allScoreData$TopScore)
  NOver50_random = nrow(filter(allScoreData, TopScore > 50))
  print('N over 50 in random sample:')
  print(NOver50_random)
  if (NOver50_random > NOver50) {
    countOver50 = countOver50 + 1
  }
  print(paste0('Iteration: ', i))
  numOver50 = c(numOver50, countOver50)
  return(allScoreData)
}

print('Proportion of runs with more genes than matched over 50: ')
print(countOver50 / 100)

amelieRandom = allRandomScores
amelieRandom$RD_ID = sapply(amelieRandom$V3, function(x)
  substr(x, start=gregexpr(pattern ='R',x)[[1]][1], stop=nchar(x)))
amelieRandom$Iteration = sapply(amelieRandom$V3, function(x)
  substr(x, start=1, stop=gregexpr(pattern ='R',x)[[1]][1]-1))
colnames(amelieRandom)[1:3] = c('Gene','AmelieScore','ID')

# Read in Amelie output
permuteData = fread('/users/nferraro/projects/rds_ase_analysis/scripts/permute.log')
kinds = which(permuteData[,2] == 'N over 50 in random sample:')
kinds2 = c(1,kinds+1)
permuteNumber = permuteData[kinds2,2]
colnames(permuteNumber) = 'Number'
permuteNumber$Permutation = 'Random'
permuteNumber$HPO_Terms = 'Random'
permuteNumber = rbind(permuteNumber, list(16,'Random', 'Matched'))
permuteNumber$Number = as.numeric(permuteNumber$Number)
rcols = c('black', 'red')
permuteNumber$HPO_Terms = factor(permuteNumber$HPO_Terms, levels=c('Random','Matched'))
ggplot(permuteNumber, aes(x=Permutation,y=Number,Group=HPO_Terms)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1, aes(fill=HPO_Terms)) +
  scale_fill_manual(values=rcols) + theme_bw() +
  ylab('Number of genes with Amelie score > 50') +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        axis.text.x=element_blank(),
        legend.title=element_text(size=20),
        legend.text=element_text(size=18))

allRDFiles = dir('/users/nferraro/projects/rds_ase_analysis/scripts/', '*rank.txt',full=T)

allRandomScoreData = do.call(rbind, lapply(allRDFiles, function(x) getAllScores(x)))

allRandomScoreData$TopScore = as.numeric(allRandomScoreData$TopScore)
allRandomScoreData$TopScore = sapply(allRandomScoreData$TopScore, function(x)
  ifelse(is.na(x), 0, x))

# Removing RD058 as a global ASE outlier (has many more outliers than other cases)
filteredRandomScoreData = filter(allRandomScoreData, TopScore != '') %>%
  mutate(HPO_terms = 'Random') %>%
  filter(!grepl('RD058',ID)) %>% group_by(ID) %>%
  mutate(MedScore = max(TopScore)) %>% 
  mutate(PropNonZero = length(which(TopScore > 0))/n()) %>%
  sample_n(1) %>%
  ungroup()

rdsCaseScores = allScoreData
rdsCaseScores$TopScore = as.numeric(rdsCaseScores$TopScore)
rdsCaseScores$TopScore = sapply(rdsCaseScores$TopScore, function(x)
  ifelse(is.na(x), 0, x))
filteredCaseScores = filter(rdsCaseScores, !is.na(TopScore)) %>%
  mutate(HPO_terms = 'Matched') %>% filter(ID != 'RD058') %>%
  filter(ID != '0RD058') %>% group_by(ID) %>%
  mutate(MedScore = max(TopScore)) %>% 
  mutate(PropNonZero = length(which(TopScore > 0))/n()) %>%
  sample_n(1) %>%
  ungroup()

colnames(filteredCaseScores)[1] = 'Gene'
allScores = rbind(filteredRandomScoreData, filteredCaseScores)
allScores$TopScore = as.numeric(allScores$TopScore)

# Wilcoxon rank sum test for difference in score distribution - alternative that Matched is greater
scorePval = wilcox.test(filter(allScores, HPO_terms == 'Matched')$PropNonZero,
                        filter(allScores, HPO_terms == 'Random')$PropNonZero,
                        alternative='g')$p.value

# Plot in figure 4 - Laure updated color/theme elsewhere
allScores$HPO_terms = factor(allScores$HPO_terms, levels=c('Random','Matched'))
scoresPlot = ggplot(allScores, aes(x=HPO_terms,y=PropNonZero)) +
  geom_violin() + geom_boxplot(width=0.1) + xlab('') +
  theme_bw() + ylab('Proportion of genes with Amelie score > 0 per sample') +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18),
        legend.title=element_text(size=20),
        legend.text=element_text(size=18)) +
  annotate("text",x=1.5,y=0.675,label=paste0('p = ',round(scorePval,3)),cex=10)
