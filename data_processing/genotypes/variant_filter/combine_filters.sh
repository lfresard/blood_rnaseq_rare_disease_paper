#!/bin/bash

variant_filter=~/projects/udn/variant_filter/vcf_processed/filter_byind/
site_filter=~/projects/udn/variant_filter/vcf_processed/udn_site_filters.txt
temp_site=~/projects/udn/variant_filter/vcf_processed/temp_site.txt

# for parallel quoting is necessary for ">" otherwise will be interpreted by bash and not feed to parallel as arguments

convert_contig=~/projects/tablejoin/contig/convert_contig.sh


bash ${convert_contig} ${site_filter} > ${temp_site}

wait

ls -1 ${variant_filter}/*.filter.txt | awk -F"\t" -v vf=${variant_filter} 'BEGIN{OFS="\t"}{match($1, /([^\/]+)\.filter\.txt/, arr); print $1,vf"/"arr[1]".filter_temp.txt"}' \
        | parallel --jobs 20 --colsep "\t" --replace {} bash ${convert_contig} {1} ">" {2}

wait

ls -1 ${variant_filter}/*.filter_temp.txt | awk -F"\t" -v vf=${variant_filter} 'BEGIN{OFS="\t"}{match($1, /([^\/]+)\.filter_temp\.txt/, arr); print $1,vf"/"arr[1]".filter_varsite.txt"}' \
	| parallel --jobs 20 --colsep "\t" --replace {} 'tbljoin -k"chr_numeric","pos" -l -n"NA"' {1} ${temp_site}  ">" {2}



rm ${temp_site}
rm ${variant_filter}/*.filter_temp.txt








