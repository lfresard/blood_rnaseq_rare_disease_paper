#!/bin/bash


## set the folder where the input and intermediate files will be saved
export RAREVARDIR=~/projects/udn/rv_annotation/output_feature

## all bed.gz, vcf.gz files are indexed by tabix, with associated .tbi files

#### [step 0]
## tabix indexed WGS vcf file
WGS_VCF_FILE=~/projects/udn/process_vcf/vcf_combined/udn_homogenized_short.vcf.gz
# /users/xli6/data/xin/udn/vcf/all_homogenized_short.vcf.gz
## VEP annotation
VEP_ANNOTATION=~/projects/udn/vep/output/udn_vep.vcf.gz
## computed 1KG allele frequency
AF_1KG_EUR=~/data_scratch/1000genomes/AF/EUR.chrALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.frq.gz
AF_1KG=~/data_scratch/1000genomes/AF/ALL.chrALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.frq.gz
## DANN annotation
DANN_ANNOTATION=~/data_scratch/gtex/river/RIVER_input/annotation/DANN_whole_genome_SNVs.tsv.bgz
## CADD annotation
CADD_SNP=~/data_shared/CADD/whole_genome_SNVs_inclAnno.tsv.gz
CADD_INDEL=~/data_shared/CADD/InDels_inclAnno.tsv.gz
## other annotations
CHROMHMM=~/data_scratch/gtex/river/RIVER_input/annotation/wgEncodeBroadHmmGm12878HMM.bed.gz
PHYLOP=~/data_scratch/phylop/phyloP100way/hg19.100way.phyloP100way.bedgraph.gz

FILTERED_VCF=${RAREVARDIR}/data/wgs/temp_filtered.vcf.gz
FILTERED_VCF=${WGS_VCF_FILE}


#### [Step 1] calculate af in study cohort itself
# this is special for UDN, as it is combined from multiple individual vcf, and ./. should be filled with 0/0 for calculating af vcftools
# then use vcftools --freq
MISFIL_VCF=~/projects/udn/process_vcf/vcf_combined/udn_homogenized_short_filmisw00.vcf.gz
AF_FILE=~/projects/udn/process_vcf/vcf_combined/udn_homogenized_short.vcf.frq.gz


IND_LIST=$(vcf-query -l ${WGS_VCF_FILE})


extract_rvsite=~/projects/rv_annotation/extract_rvsites_ByInd.py
extract_score=~/projects/rv_annotation/extract_scores_combined.py


# vcf already filtered this step skipped
#### [step 2]
## retain high quality variant calls (VQSLOD = PASS) were considered, only sites having <= 10 individuals in terms of missing genotypes were considered, and only autosomes were considered, and only European subjects were considered.
<<skipped
vcftools --gzvcf ${WGS_VCF_FILE} --keep ${IND_LIST} --remove-filtered-all --max-missing-count 10 --not-chr X --recode --recode-INFO-all --stdout | bgzip -c > ${FILTERED_VCF}
tabix -p vcf ${FILTERED_VCF}
skipped


#### [Step 3] calculate af within the study set
<<skipped
vcftools --gzvcf ${MISFIL_VCF} --keep <(echo ${IND_LIST} | tr ' ' '\n') --freq --stdout | bgzip -c > ${AF_FILE}
tabix -p vcf ${AF_FILE}
skipped


#### [Step 4] gerate the regions around gene, where the variants will be extracted
GENE_REGION=~/data_scratch/gtex/river/RIVER_input/preprocessing/rvsite/region.tss10k.txt

#### [Step 5] In each subject of interest, extract a list of individual-specific rare variant sites based on AFs of both GTEx and EUR 1K population
<<skipped
count_ind=0
for ID in ${IND_LIST}
do
count_ind=$(( $count_ind + 1 ))
echo "\
cat ${GENE_REGION} | python ${extract_rvsite} -n $count_ind --id $ID --WGSvcf_in ${FILTERED_VCF} --AFvcf_in ${AF_FILE} --AFpop_in ${AF_1KG} --AFcutoff 0.02 --site_out ${RAREVARDIR}/data/indiv/${ID}.${count_ind}.rvsites.txt"
done
skipped


#### [Step 6] Extract all the features simulataneously (CADD, chromHMM, phylop, DANN).
count_ind=0
for ID in ${IND_LIST}
do
count_ind=$(( $count_ind + 1))
echo "\
cat ${RAREVARDIR}/data/indiv/${ID}.$count_ind.rvsites.txt | python ${extract_score} -n $count_ind --id $ID --af_in ${AF_FILE} --wgs_in ${FILTERED_VCF} --anno_in ${VEP_ANNOTATION} --cadd_in ${CADD_SNP} --cadd_indel_in ${CADD_INDEL}  --dann_in ${DANN_ANNOTATION} --chromHMM_in ${CHROMHMM} --phylop_in ${PHYLOP} --score_out ${RAREVARDIR}/data/score/${ID}.${count_ind}.score.nuc.txt"
done

# combine results
<<skipped
realpath ${RAREVARDIR}/data/score/*.score.nuc.txt | xargs tblcat > udn.score.nuc.txt
skipped

