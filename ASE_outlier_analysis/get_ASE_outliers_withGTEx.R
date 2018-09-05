library(data.table)
library(dplyr)

data_dir = '/users/nferraro/data/rds_data/'

# Read in the output from collapse_GTEx_ASE.R
ase.all = fread(paste0(data_dir, 'gtexV7_ase_filteredSitesForRDS_maxPerTissue.txt'), data.table=F)
ase.all = ase.all[,c(1,5,6,4,7)]
colnames(ase.all) = c('SampID', 'TotalReads', 'RefRatio', 'GeneName', 'CPOS')
rds.ase.data = fread(paste0(data_dir, 'ensembleids.ase.filtered.txt'))

# Check that all sites are found in GTEx
rds.ase.data = filter(rds.ase.data, RefRatio != 0 & RefRatio != 1) %>%
 mutate(CPOS=paste(Chr,POS,sep='_')) %>% filter(CPOS %in% ase.all$CPOS) %>%
  select(SampID,Gene,TotalReads,RefRatio,CPOS)

rds.ase.data$GeneName = sapply(rds.ase.data$Gene, function(x)
  strsplit(x, '[.]')[[1]][1])

rds.ase.data = rds.ase.data %>% filter(GeneName %in% ase.all$GeneName)
ase.all = ase.all %>% filter(CPOS %in% rds.ase.data$CPOS) %>%
  filter(GeneName %in% rds.ase.data$GeneName)

rds.ase.data = rds.ase.data[,c(1,6,3,4,5)]
ase.all = ase.all[,c(1,4,2,3,5)]
# Scale the reference ratios for all sites within each gene to obtain z-scores
rds_gtex_single_ase = rbind(ase.all, rds.ase.data) %>% group_by(GeneName) %>%
  mutate(RatioScaled = scale(RefRatio)) %>% ungroup()

write.table(rds_gtex_single_ase, file=paste0(data_dir, 'RDS_GTEX_ASE_zscores.txt'), sep='\t', quote=F, row.names=F)

rds_gtex_single_outliers = rds_gtex_single_ase %>% group_by(GeneName) %>%
  mutate(UpRank = rank(-RatioScaled)) %>% mutate(DownRank = rank(RatioScaled)) %>%
  mutate(BothRank = rank(-abs(RatioScaled))) %>% ungroup()

### Look at number of outliers per individual
mostExtreme = filter(rds_gtex_single_outliers, BothRank == 1)
colnames(mostExtreme)[1] = 'ID'
inds2 = as.data.frame(table(mostExtreme$ID)) %>%
  mutate(Source = ifelse(grepl('RD', Var1), 'RDS', 'GTEX')) %>%
  mutate(Status = ifelse(Var1 %in% conids, 'Case', 'Control')) %>%
  mutate(Set = paste(Source, Status,sep='_')) %>%
  mutate(Threshold = 2)
inds3 = as.data.frame(table(filter(mostExtreme, abs(RatioScaled) > 3)$ID)) %>%
  mutate(Source = ifelse(grepl('RD', Var1), 'RDS', 'GTEX')) %>%
  mutate(Status = ifelse(Var1 %in% conids, 'Case', 'Control')) %>%
  mutate(Set = paste(Source, Status,sep='_')) %>%
  mutate(Threshold = 3)
inds4 = as.data.frame(table(filter(mostExtreme, abs(RatioScaled) > 4)$ID)) %>%
  mutate(Source = ifelse(grepl('RD', Var1), 'RDS', 'GTEX')) %>%
  mutate(Status = ifelse(Var1 %in% conids, 'Case', 'Control')) %>%
  mutate(Set = paste(Source, Status,sep='_')) %>%
  mutate(Threshold = 4)
inds5 = as.data.frame(table(filter(mostExtreme, abs(RatioScaled) > 5)$ID)) %>%
  mutate(Source = ifelse(grepl('RD', Var1), 'RDS', 'GTEX')) %>%
  mutate(Status = ifelse(Var1 %in% conids, 'Case', 'Control')) %>%
  mutate(Set = paste(Source, Status,sep='_')) %>%
  mutate(Threshold = 5)

allCounts = rbind(inds2,inds3,inds4,inds5)
allCounts$Threshold = factor(allCounts$Threshold, levels=c(2,3,4,5))
write.table(allCounts, file=paste0(data_dir, 'RDS_ASE_outlier_counts.txt'), sep='\t', quote=F, row.names=F)




