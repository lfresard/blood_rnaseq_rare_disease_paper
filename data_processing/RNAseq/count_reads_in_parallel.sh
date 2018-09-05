#!/bin/bash
# This script counts the number of reads we have sequenced for this project

find /srv/scratch/restricted/rare_diseases/data/fastq/ \
    | grep fastq.gz \
    | sort \
    | parallel  "./count_reads.sh {}" \
    > /srv/scratch/restricted/rare_diseases/analysis/fastq_analysis/raw_read_counts.txt

