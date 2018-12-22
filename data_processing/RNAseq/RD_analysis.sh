#!/bin/bash

batch_number=$1

#run quantification piepeline for each sample


FASTQ_DIR=<fastq_dir>/batch${batch_number}
cutadapt_script=trim_adapters_150bpreads.sh #include path to script
STAR_script=STAR_alignment_rare_disease.sh #include path to script

BAM_DIR=<bam_dir>/batch${batch_number}
FILTER_script=Filters_bam_uniq_mq30_duplicates_rare_disease_transcriptome.sh #include path to script
EXP_DIR=<quantification_dir>rsem
RSEM_script=rsem_calculate_expression_raredisease.sh  #include path to script

date



echo 'start trimming'
cd ${FASTQ_DIR}

ls *merge_R1.fastq.gz | sed 's/_/\t/'| awk '{print $1}' |awk -v fastq_dir=$FASTQ_DIR 'BEGIN{OFS="\t"}{print fastq_dir"/"$1"_merge_R1.fastq.gz",fastq_dir"/"$1"_merge_R1.trimmed.fastq.gz", fastq_dir"/"$1"_merge_R2.fastq.gz",fastq_dir"/"$1"_merge_R2.trimmed.fastq.gz", fastq_dir"/"$1"_merge_short_reads.fastq.gz",fastq_dir"/"$1"_fastq_trimadapters.log"}' | \
parallel --jobs 15 --col-sep "\t" "${cutadapt_script} {1} {2} {3} {4} {5} {6}" 
echo 'done trimming'


echo "start mapping"
echo ""


ls *merge_R1.trimmed.fastq.gz | sed 's/_/\t/'| awk '{print $1}' |awk -v fastq_dir=$FASTQ_DIR -v bam_dir=$BAM_DIR 'BEGIN{OFS="\t"}{print fastq_dir"/"$1"_merge_R1.trimmed.fastq.gz", fastq_dir"/"$1"_merge_R2.trimmed.fastq.gz",bam_dir"/"$1}' |\
	parallel --jobs 5 --col-sep "\t" "${STAR_script} {1} {2} {3} $<index_dir>/STAR_INDEX_OVERHANG_150/ 150"

wait

STAR --genomeDir $<index_dir>/STAR_INDEX_OVERHANG_150 --genomeLoad Remove

echo "mapping finished"
echo ""


echo "remove fastq files"
echo ""
#rm $FASTQ_DIR/*merge_R1.trimmed.fastq.gz
#rm $FASTQ_DIR/*merge_R2.trimmed.fastq.gz


echo "start Bam filter for mapping quality and PCR duplicates"
echo ""
cd $BAM_DIR


ls ${BAM_DIR}/*.Aligned.toTranscriptome.out.bam | parallel --jobs 10 --col-sep "\t" "${FILTER_script} {1}"

echo "Filter done"
echo ""




echo "start quantifying batch ${batch_number}"
cd $BAM_DIR

ls *.Aligned.toTranscriptome.out_mapq30_sorted_dedup_byread.bam | sed 's/\./\t/'|awk '{print $1}' |awk -v bam_dir=$BAM_DIR  'BEGIN{OFS="\t"}{print bam_dir"/"$1".Aligned.toTranscriptome.out_mapq30_sorted_dedup_byread.bam",$1}' |   \
	parallel   --jobs 10 --col-sep "\t" "bash ${RSEM_script} {1} {2}"


echo "Remove filtered bam files"
echo ""

rm *.Aligned.out_mapq30_sorted_dedup.ba*
rm *.Aligned.toTranscriptome.out_mapq30_sorted_dedup_byread.ba*


echo "finished quantifying batch ${batch_number}"
echo ""
date

