#!/bin/R.3.3.2
#LF

# This script take the metadata file and the tissue of interest as an input
# It creates a directory for each case
# in each of those directory, it creates a list of case and controls


###
### Example usage: /users/lfresard/R_3.3.2/bin/Rscript make_list_case_controls.R --meta <metadata_file> --tissue <tissue_to_analyze> --outdir <path_to_results> --juncdir <path_to_junc_files> --freeze FALSE
###





#--- LIBRARIES

library(optparse)
library(dplyr)
library(tximport)
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
library(preprocessCore)
library(biobroom)


#--- OPTION PARSER
option_list = list(
  make_option(c("-m", "--meta"), type="character", default="/srv/scratch/restricted/rare_diseases/data/metadata/metadata.tsv", help="metadata file name", metavar="character"),
  make_option(c("-t", '--tissue'), type="character", default=NULL, help="tissue to analyze", metavar="character"),
  make_option(c("-o", '--outdir'), type="character", default=".", help="path to output [default= %default]", metavar="character"),
  make_option(c("-j", '--juncdir'), type="character", default=".", help="path to junctions [default= %default]", metavar="character"),
  make_option(c("-f", "--freeze"), type="logical", default=FALSE, help="wether the analysis is conducted on the freeze", metavar="character"),
  make_option(c("-d", "--DGN"), type="logical", default=FALSE, help="wether the analysis is conducted on with or without DGN", metavar="character"),
  make_option(c("-p", "--PIVUS"), type="logical", default=FALSE, help="wether the analysis is conducted on with or without PIVUS", metavar="character")
) 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#--- VARIABLES

metadata=opt$meta
tissue=opt$tissue
output_dir=opt$outdir
junction_dir=opt$juncdir
freeze_analysis=opt$freeze
DGN_included=opt$DGN
PIVUS_included=opt$PIVUS
sampling_mode=opt$sampling_mode
#sampling_nb=opt$sampling_number


# Filter metadata depending if studying freeze data or not.
# If positive, integrate DGN samples in analysis
if(DGN_included==TRUE){
	DGN_samples=data.frame(sample_id=gsub(".SJ.out_filtered_uniq_10.tab", "", list.files('<junc_dir>/all_filteredjunc', pattern="LD")),affected_status=rep("Control", length(list.files('<junc_dir>/all_filteredjunc', pattern="LD"))))
} else {DGN_samples=NULL}

if(PIVUS_included==TRUE){
	PIVUS_samples=data.frame(sample_id=gsub(".SJ.out_filtered_uniq_10.tab", "", list.files('<junc_dir>/all_filteredjunc', pattern="PIVUS")),affected_status=rep("Control", length(list.files('<junc_dir>/all_filteredjunc', pattern="PIVUS"))))
} else {PIVUS_samples=NULL}

if (freeze_analysis==TRUE){
	metadata = read_tsv(metadata) %>% filter(in_freeze == 'yes') %>% filter(source_of_RNA=="Blood")
	metadata=metadata %>% select(sample_id,affected_status)
	metadata=rbind(metadata,DGN_samples,PIVUS_samples)
} else{
	metadata = read_tsv(metadata) %>% filter(source_of_RNA == tissue & sequencing_status == 'PASSED' & is_RD=='yes')
	metadata=metadata %>% select(sample_id,affected_status)
	metadata=rbind(metadata,DGN_samples,PIVUS_samples)
}


#--- MAIN

setwd(output_dir)
junc_files=metadata %>% select(sample_id) %>% unlist %>% paste0(".SJ.out_filtered_uniq_10.tab")
to_print=as.data.frame(paste0(junction_dir,'/', junc_files))


if(DGN_included==TRUE){
	meta_filename='sample_affected_status_freeze_RD_DGN.tsv'
	filename='list_junctions_freeze_RD_DGN.txt'
} else if (PIVUS_included==TRUE){
		meta_filename='sample_affected_status_freeze_RD_PIVUS.tsv'
		filename='list_junctions_freeze_RD_PIVUS.txt'
}else {
	meta_filename='sample_affected_status_freeze_RD.tsv'
	filename='list_junctions_freeze_RD.txt'
}


	#filename_control = 'list_controls.txt'
write_tsv(metadata,meta_filename, col_names=F )
write_tsv(to_print, filename, col_names=F)





