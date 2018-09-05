#!/bin/bash
chr=$1
njob=15
gnomAD_dir=<gnomad_data_dir>
cadd_dir=<CADD_data_dir>
# for each chromosome:
	#parse the vcf files
	# get the AF from gnomAD

# At the end concatenate all chromosome back together
export chr 
ls RD*_homogenized.vcf.gz| sed 's/\_/\t/'| awk '{print $1}' | parallel --jobs $njob  \
	"vcftools --gzvcf {}_homogenized.vcf.gz --chr ${chr} --recode --stdout | bgzip -c > {}_homogenized_chr${chr}.vcf.gz"

wait

echo '[[annotation]]'> conf.toml 
echo 'file="'${gnomAD_dir}/gnomad.genomes.r2.0.2.sites.chr${chr}.vcf.gz'"' >> conf.toml 
echo 'fields = ["AF"]'>> conf.toml
echo 'names=["gnomAD_AF"]' >> conf.toml 
echo 'ops=["first"]' >> conf.toml 
echo "">> conf.toml 
echo "">> conf.toml 
echo '[[annotation]]'>> conf.toml 
echo 'file="'$cadd_dir/CADD_SNP.INDELS.sorted.vcf.gz'"'>> conf.toml 
echo 'names=["cadd_phred", "cadd_raw"]'>> conf.toml 
echo 'ops=["mean", "mean"]'>> conf.toml 
echo 'fields=["phred", "raw"]'>> conf.toml 


ls RD*_homogenized_chr${chr}.vcf.gz | sed 's/\_/\t/'| awk '{print $1}' | parallel --jobs $njob \
	"bash vcfanno.sh {} ${chr}"

wait

rm conf.toml
rm RD*_homogenized_chr${chr}.vcf.gz

echo "Done with chromosome ${chr}"
