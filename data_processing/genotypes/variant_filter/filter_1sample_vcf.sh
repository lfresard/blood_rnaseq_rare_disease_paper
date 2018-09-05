#!/bin/bash

site_filter=$1
vcf_file=$2
filtered_file=$3

temp_file=${site_filter}"_temp"

tblmap -k"chr,pos,id,ref,alt,chr_numeric,vflag,chr.1,f_miss,p_hwe" ${site_filter} | \
	awk -F"\t" 'BEGIN{OFS="\t"}{ \
		fail=$7=="FAIL" || $9 > 0.2 || $10 < 1e-6; \
		if(fail){print $1,$2} \
		}' > ${temp_file}

vcftools --gzvcf ${vcf_file} --exclude-positions ${temp_file} --recode --stdout --out ${temp_file} | vcf-sort | bgzip > ${filtered_file} 

rm ${temp_file}


