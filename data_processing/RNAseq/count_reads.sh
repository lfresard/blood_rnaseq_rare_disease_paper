#!/bin/bash
# This script simply counts the number of reads in a gzipped FASTQ file

fastq=$1

linecount=$(zcat $fastq | wc -l)
readcount=$((linecount/4))

echo -e $fastq"\t"$linecount"\t"$readcount
