sample=$1

ls ${sample}_homogenized_chr*_gnomad_cadd.vcf.gz > ${sample}_list.txt

wait


vcf-concat --files ${sample}_list.txt |bgzip -c >${sample}_homogenized_gnomad_cadd.vcf.gz
tabix ${sample}_homogenized_gnomad_cadd.vcf.gz

