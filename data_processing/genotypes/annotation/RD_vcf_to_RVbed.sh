#!/bin/bash
#LF
#Jan 2018

# This script takes annotated sample vcf with gnomAD allele frequency and CADD scores 
# filters for Rare variants
# transforms into bed file
# get genes name overlapping with rare variant

vcf_file=$1
fname=`basename $vcf_file`
sample=${fname%%_*}
bed_dir=/srv/scratch/restricted/rare_diseases/data/bed/for_freeze_exac_filt
out_file=${bed_dir}/${fname%.vcf.gz}_RV.bed
gene_bed=/srv/scratch/restricted/rare_diseases/data/annot/gencode.v19.annotation.genes.bed.gz
out_file_final=${bed_dir}/${fname%.vcf.gz}_RV_withgene.bed.gz

date
echo ""
echo "Filtering ${vcf_file} for rare variants..."
echo ""

## Filter file for rare variant
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%gnomAD_AF\t%cadd_raw\t%cadd_phred\n' $vcf_file | \
	awk '$5<=0.01 {OFS="\t"; print "chr"$1, $2-1, $2, $3, $4, $5, $6, $7}' |sort -k1,1 -k2,2n > $out_file


wait

echo "Done filtering ${vcf_file}."
echo ""

echo "Gzipping bed file for ${sample}..."
echo ""

gzip $out_file

echo "Done gzipping ${sample}."
echo ""

echo "Find gene location ${sample}"
echo ""

bedtools intersect -sorted -a ${out_file}.gz -b ${gene_bed} -wa -wb |gzip - > ${out_file_final}


echo "Done with intersection with gene bed."
echo ""

echo "Removing intermediate files."
echo ""
rm ${out_file}.gz

echo "Done with $sample"
echo ""

date
