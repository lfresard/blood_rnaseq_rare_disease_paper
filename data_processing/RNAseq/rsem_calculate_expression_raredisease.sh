#!/bin/bash
#LF
#RSEM quantification on filtered bam files (transcriptome ref) for rare disease samples
INPUT=$1 #.Aligned.toTranscriptome.out_mapq30_sorted_dedup.bam format
REF=/srv/scratch/restricted/rare_diseases/data/quantification/rsem/rsem_reference/hg19_RSEM
WORK_DIR=/srv/scratch/restricted/rare_diseases/data/quantification/rsem
N_THREADS=4
SAMPLE_NAME=${WORK_DIR}/$2
log_file=${SAMPLE_NAME}_rsem.log


date>$log_file


cmd="/users/lfresard/bin/RSEM-1.2.21/rsem-calculate-expression \
	-p $N_THREADS \
	--paired-end \
	--no-bam-output \
	--forward-prob 0.5 \
	--seed 12345 \
	--bam $INPUT\
	$REF \
	$SAMPLE_NAME"

echo $cmd >>$log_file
eval $cmd >>$log_file 2>&1

date>> $log_file

rsem_command="$rsem \
        -p $nthreads \
        --paired-end \
        --no-bam-output \
        --forward-prob 0.5 \
        --seed $seed \
        --bam ${sample}.resorted.bam \
        $rsem_reference \
        $sample"
