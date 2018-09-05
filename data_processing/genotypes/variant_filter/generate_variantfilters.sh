#!/bin/bash


# format varies across sites, samples and institutions

# bcftools does not work on UDN vcfs, format not correct

parse_vcf_cmd=~/projects/udn/variant_filter/parse_1sample_vcf.sh
output_filter=~/projects/udn/variant_filter/vcf_processed/filter_byind/
convert_contig=~/projects/tablejoin/contig/convert_contig.sh

instns=("CGS" "CHEO" "UDN")


date > ${output_filter}/udn_vcffield_FORMAT.txt
date > ${output_filter}/udn_vcffield_FILTER.txt
date > ${output_filter}/udn_contig.txt

for ins in "CGS" "CHEO" "UDN"
do
	for i in $(ls -1 vcf_all/${ins}/RD*.vcf.gz)
	do
		zcat $i | awk -F"\t" -vins=${ins} -vi=$i 'BEGIN{OFS="\t"}$2 ~ /[0-9]+/{print ins,i,$9;}' | sort | uniq >> ${output_filter}/udn_vcffield_FORMAT.txt &
		zcat $i | awk -F"\t" -vins=${ins} -vi=$i 'BEGIN{OFS="\t"}$2 ~ /[0-9]+/{print ins,i,$7;}' | sort | uniq >> ${output_filter}/udn_vcffield_FILTER.txt &
		zcat $i | awk -F"\t" -vins=${ins} -vi=$i 'BEGIN{OFS="\t"}$2 ~ /[0-9]+/{print ins,i,$1;}' | sort | uniq >> ${output_filter}/udn_contig.txt &
		prefix=$(basename ${i});
		echo "bash ${parse_vcf_cmd} $i > ${output_filter}/${prefix}.filter.txt"
	done
done









