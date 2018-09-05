#!/bin/bash


# format varies across sites, samples and institutions

# bcftools does not work on UDN vcfs, format not correct

parse_vcf_cmd=~/projects/udn/variant_filter/parse_1sample_vcf.sh
filter_vcf_cmd=~/projects/udn/variant_filter/filter_1sample_vcf.sh
output_filter=~/projects/udn/variant_filter/vcf_processed/filter_byind/
convert_contig=~/projects/tablejoin/contig/convert_contig.sh

instns=("CGS" "CHEO" "UDN")



for ins in "CGS" "CHEO" "UDN"
do
	for i in $(ls -1 vcf_all/${ins}/RD*.vcf.gz)
	do
		filename=$(basename ${i});	
		prefix=$(basename ${i} .vcf.gz);
		echo "bash ${filter_vcf_cmd}  ${output_filter}/${filename}.filter_varsite.txt $i ${output_filter}/${prefix}.filtered.vcf.gz"
	done
done









