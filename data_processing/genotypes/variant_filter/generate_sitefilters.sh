#!/bin/bash


vcf_combined=~/projects/udn/variant_filter/vcf_processed


vcftools --gzvcf ${vcf_combined}/all_homogenized_short.vcf.gz --missing-site --stdout | bgzip > ${vcf_combined}/udn_vcf_missing.txt.gz
vcftools --gzvcf ${vcf_combined}/all_homogenized_short.vcf.gz --hardy --stdout | bgzip > ${vcf_combined}/udn_vcf_hwe.txt.gz


<<skipped
# another option is to use bedtools to link
# -f 1 to force fully overlapping
bedtools intersect -a <(vcftools --gzvcf vcf_combined/all_homogenized_short.vcf.gz --missing-site --stdout | awk 'BEGIN{OFS="\t"}{print $1,$2,$2,$6}' | head) \
	-b <(vcftools --gzvcf vcf_combined/all_homogenized_short.vcf.gz --hardy --stdout | awk 'BEGIN{OFS="\t"}{print $1,$2,$2,$6}' | head) \
	-wa -wb -loj -f 1
skipped

# tbljoin converts all to lower case column names
tbljoin -k"chr","pos" \
	-lr \
	${vcf_combined}/udn_vcf_missing.txt.gz \
	${vcf_combined}/udn_vcf_hwe.txt.gz \
	| awk 'BEGIN{OFS="\t"}{print $1,$2,$6,$10;}' > ${vcf_combined}/udn_site_filters.txt


rm ${vcf_combined}/udn_vcf_missing.txt.gz
rm ${vcf_combined}/udn_vcf_hwe.txt.gz




