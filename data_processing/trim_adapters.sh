#!/bin/bash
# This script trim adapters from fastqs (75bp reads)

# Set variables
IN1=$1
OUT1=$2
IN2=$3
OUT2=$4
Too_short=$5
log_file=$6

###
### TrueSeq Sequences from Kevin
###
INDEXED_ADAPTER_PREFIX=GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
UNIVERSAL_ADAPTER=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT 



###
### Paired-end trimming of adapters
###
date >$log_file
echo "trimming adapters" >>$log_file
echo "INDEXED ADAPTER PREFIX:">>$log_file
echo $INDEXED_ADAPTER_PREFIX>>$log_file

echo "UNIVERSAL ADAPTER:">>$log_file
echo $UNIVERSAL_ADAPTER >> $log_file

cmd="cutadapt \
        -a $INDEXED_ADAPTER_PREFIX \
        -A $UNIVERSAL_ADAPTER \
        -o $OUT1 -p $OUT2 \
        -m 20 \
        --too-short-paired-output=$Too_short \
        $IN1 $IN2 "

echo $cmd >>$log_file
eval $cmd >> $log_file 2>&1
echo "Trimming done" >>$log_file

date >>$log_file








