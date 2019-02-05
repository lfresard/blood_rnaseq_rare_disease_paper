cd $WORK_DIR

# 1. Get sample IDs from outlier file names and write to a text file sample_id.txt

# 2. Sort expression outlier BED files
parallel --verbose --xapply -a sample_id.txt \
  'sort -k1,1 -k2,2n outliers/{1}_outliers.bed > outliers/{1}_outliers_sorted.bed'

# 3. Filter VCF for variants with MAF < 0.01, then bedtools closest command to get variants matched to closest expression outlier genes
parallel --verbose --xapply -a sample_id.txt \
  'bcftools query -f '\''%CHROM\t%POS\t%REF\t%ALT\t%gnomAD_AF\t%cadd_raw\t%cadd_phred\n'\'' <path_to_vcf>/{1}_homogenized_gnomad_cadd.vcf.gz | \
    awk '\''$5<=0.01 {OFS="\t"; print $1, $2-1, $2, $3, $4, $5, $6, $7}'\'' | sort -k1,1 -k2,2n | \
    bedtools window -w 10000 -a outliers/{1}_outliers_sorted.bed -b stdin | \
    awk -F "\t" '\''{OFS="\t"; print $5, $6, $4, $8, $12, $13, $14}'\'' > outliers/{1}_tmp.txt'

# 4. Concatenate files in to one
cd outliers
ls *_tmp.txt > ../output_files.txt
cd ..

while read FILE_NAME
do
   cat outliers/${FILE_NAME} >> outliers/outliers_rare_var_combined_10kb.txt
done < output_files.txt


