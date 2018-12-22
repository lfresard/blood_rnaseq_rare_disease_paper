#!/bin/bash
for i in `seq 1 22`;
do
    bash vcf_anno_chr.sh $i
done    
       
wait

bash vcf_anno_chr.sh X

wait

ls *_gnomad.vcf.gz| sed 's/\_/\t/'| awk '{print $1}' | parallel --jobs $njob  \
	
ls RD*_homogenized_chr*_gnomad_cadd.vcf.gz | sed 's/\_/\t/'| awk '{print $1}' | parallel --jobs 3 \
	"bash concatenate_gnomadres.sh {}"


wait

rm RD*_homogenized_chr*_gnomad_cadd.vcf.gz
