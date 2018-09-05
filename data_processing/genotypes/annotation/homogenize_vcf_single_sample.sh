#!/bin/bash
# June 2018
# Laure Fresard
# Transform vcfs filtered by Xin into homogeneous format for further analysis

date


vcf=$1 # /srv/scratch/restricted/rare_diseases/data/vcfs/before_merge
sample=$2 # name to put in the resulting VCF column
tempdir=$3 #"/srv/scratch/restricted/rare_diseases/data/vcfs/intermediates"

# output sorted vcf body
header="${tempdir}/${sample}_header.txt"
body="${tempdir}/${sample}_body.txt"
body2="${tempdir}/${sample}_body2.txt" #samplify info field before merge


# create homogenized files
echo "create homogenized files"
echo ""


# print out the header
timestamp=$(date +"%m-%d-%Y")
echo "##fileformat=VCFv4.1" > $header
echo "##fileDate=$timestamp" >> $header
echo "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">" >> $header
echo "##FORMAT=<ID=GT,Number=1,Type=Integer,Description=\"Genotype\">" >> $header
echo "##contig=<ID=1,length=249250621,assembly=GRCh37>" >> $header
echo "##contig=<ID=2,length=243199373,assembly=GRCh37>" >> $header
echo "##contig=<ID=3,length=198022430,assembly=GRCh37>" >> $header
echo "##contig=<ID=4,length=191154276,assembly=GRCh37>" >> $header
echo "##contig=<ID=5,length=180915260,assembly=GRCh37>" >> $header
echo "##contig=<ID=6,length=171115067,assembly=GRCh37>" >> $header
echo "##contig=<ID=7,length=159138663,assembly=GRCh37>" >> $header
echo "##contig=<ID=8,length=146364022,assembly=GRCh37>" >> $header
echo "##contig=<ID=9,length=141213431,assembly=GRCh37>" >> $header
echo "##contig=<ID=10,length=135534747,assembly=GRCh37>" >> $header
echo "##contig=<ID=11,length=135006516,assembly=GRCh37>" >> $header
echo "##contig=<ID=12,length=133851895,assembly=GRCh37>" >> $header
echo "##contig=<ID=13,length=115169878,assembly=GRCh37>" >> $header
echo "##contig=<ID=14,length=107349540,assembly=GRCh37>" >> $header
echo "##contig=<ID=15,length=102531392,assembly=GRCh37>" >> $header
echo "##contig=<ID=16,length=90354753,assembly=GRCh37>" >> $header
echo "##contig=<ID=17,length=81195210,assembly=GRCh37>" >> $header
echo "##contig=<ID=18,length=78077248,assembly=GRCh37>" >> $header
echo "##contig=<ID=19,length=59128983,assembly=GRCh37>" >> $header
echo "##contig=<ID=20,length=63025520,assembly=GRCh37>" >> $header
echo "##contig=<ID=21,length=48129895,assembly=GRCh37>" >> $header
echo "##contig=<ID=22,length=51304566,assembly=GRCh37>" >> $header
echo "##contig=<ID=X,length=155270560,assembly=GRCh37>" >> $header
echo "##contig=<ID=Y,length=59373566,assembly=GRCh37>" >> $header
echo "##contig=<ID=MT,length=16571,assembly=GRCh37>" >> $header
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sample" >> $header


# output sorted vcf body
zcat $vcf | grep -v "^#" | sed s/chr// | sed s/M/MT/ | awk 'BEGIN{FS="\t";OFS="\t"} length($1) < 3 {print $1, $2, ".", $4, $5, ".", $7, ".",$9, $10}' > $body
zcat $vcf | grep -v "^#" | sed s/chr// | sed s/M/MT/ | awk 'BEGIN{FS="\t";OFS="\t"} length($1) < 3 {print $1, $2, ".", $4, $5, ".", "PASS", "NS=1","GT", substr($10,1,3)}' | grep -v 1/2 |sort -k1,1d -k2,2n >$body2

wait

# combine final vcf header and body
finalvcf="${tempdir}/${sample}_homogenized.vcf"
finalvcf2="${tempdir}/${sample}_homogenized_short.vcf"

cat $header $body > $finalvcf
cat $header $body2> $finalvcf2

wait
echo "bgzip and tabix homogenized vcf"
echo ""


wait
# bgzip and tabix final vcf
bgzip -f $finalvcf
tabix $finalvcf.gz

bgzip -f $finalvcf2
tabix $finalvcf2.gz


# remove intermediate files
echo "clean intermediate files"
echo ""

rm $header $body $body2 



echo "Done"
echo ""
echo ""
date

