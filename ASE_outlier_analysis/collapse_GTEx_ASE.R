## Looking into matching of ASE with outlier expression tissue by tissue.
## Testing the extent to which tissues with extreme expression also show ASE enrichments.

library(data.table)
library(dplyr)
library(ggplot2)
require(RColorBrewer)
library(argparse)

parser = ArgumentParser()
parser$add_argument('--asedir', help = 'Prefix of outlier files.')
parser$add_argument('--datadir', help = 'Window name.')
parser$add_argument('--rdsdir', default = 0.01, help = 'FDR threshold.')
args = parser$parse_args()

ASEDIR = args$asedir
data_dir = args$datadir
rds_dir = args$rdsdir

########### FUNCTIONS - adopted from GTEx ASE analysis (Emily Tsang)
## Function to read in ASE data from a file and apply certain QC filters
## Input: file name
## Output: data.table of filtered ASE data
read.ase = function(filename) {
    ase.data = fread(paste('zcat', filename))
    nall = nrow(ase.data)
    # Quality control filter + remove absolutes
    ase.data = filter(ase.data,
                      LOW_MAPABILITY == 0 & MAPPING_BIAS_SIM == 0 & GENOTYPE_WARNING == 0 &
                      TOTAL_COUNT >=20, REF_RATIO != 0, REF_RATIO != 1)
    # note the data don't include sex chromosomes, so not including that filter
    nfilt = nrow(ase.data)
    cat(filename, '\n')
    cat('number of rows before filters:', nall, '\n')
    cat('number of rows after filters:', nfilt, '\n')
    cat('percentage of rows after filters:', nfilt/nall, '\n')
    return(ase.data)
}


## Read in all the v7 ASE data
ase.files = list.files(ASEDIR, full.names = TRUE)

ase.list = lapply(ase.files, read.ase)
ase.all = rbindlist(ase.list)
rm(ase.list)

ase.all = ase.all %>% select(CHR,POS,VARIANT_ID,SUBJECT_ID,TISSUE_ID,
                             REF_COUNT,TOTAL_COUNT,REF_RATIO,GENE_ID) %>%
  mutate(CPOS = paste(paste0('chr', CHR),POS,sep='_')) %>% group_by(SUBJECT_ID) %>%
  mutate(NS=length(unique(CPOS))) %>% ungroup() %>% mutate(Source = 'GTEX') %>% mutate(Status='Control') %>%
  select(SUBJECT_ID,CHR,POS,GENE_ID,TOTAL_COUNT,REF_RATIO,CPOS,NS,Status,Source)

rds.ase.data = fread(paste0(data_dir, 'ensembleids.ase.filtered.txt'))

metadata = fread(paste0(rds_dir, 'metadata/2018_06_12_Rare_Disease_Metadata.tsv'), data.table=F)
rownames(metadata) = metadata$sample_id
freeze_metadata = filter(metadata, in_freeze == 'yes')
cids = filter(freeze_metadata, affected_status == 'Case')$sample_id
conids = filter(freeze_metadata, affected_status == 'Control')$sample_id

rds.ase.data = filter(rds.ase.data, RefRatio != 0 & RefRatio != 1) %>%
 mutate(CPOS=paste(Chr,POS,sep='_')) %>% filter(CPOS %in% ase.all$CPOS)

# Take most extreme ASE per tissue per gene-individual
ase.all = ase.all %>% filter(CPOS %in% rds.ase.data$CPOS) %>% mutate(AbsRatio = abs(0.5 - REF_RATIO)) %>%
  group_by(CPOS,SUBJECT_ID) %>%
  top_n(1,AbsRatio) %>% sample_n(1) %>% ungroup()

write.table(ase.all, file=paste0(data_dir, 'gtexV7_ase_filteredSitesForRDS_maxPerTissue.txt'), sep='\t', quote=F, row.names=F)
