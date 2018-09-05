#!/bin/bash
#January 2018
#LF



#run quantification piepeline for each sample of PIVUS project


FASTQ_DIR=/srv/scratch/restricted/rare_diseases/data/pivus/fastq
STAR_script=/users/lfresard/repos/rare_disease/scripts/mapping/STAR_alignment_rare_disease.sh #include path

BAM_DIR=/srv/scratch/restricted/rare_diseases/data/pivus/bams
FILTER_script=/users/lfresard/repos/rare_disease/scripts/mapping/Filters_bam_uniq_mq30_duplicates_rare_disease.sh
SORT_SCRIPT=/users/lfresard/repos/rare_disease/scripts/mapping/sort_bam_readname.sh

EXP_DIR=/srv/scratch/restricted/rare_diseases/data/quantification/rsem
RSEM_script=/users/lfresard/repos/rare_disease/scripts/quantification/rsem_calculate_expression_raredisease_pivus.sh #with PATH

date
echo "start mapping"
echo ""
cd ${FASTQ_DIR}

ls *merge_R1.trimmed.fastq.gz | sed 's/_/\t/'| awk '{print $1}' |awk 'BEGIN{OFS="\t"}{print "/srv/scratch/restricted/rare_diseases/data/pivus/fastq/"$1"_merge_R1.trimmed.fastq.gz", "/srv/scratch/restricted/rare_diseases/data/pivus/fastq/"$1"_merge_R2.trimmed.fastq.gz","/srv/scratch/restricted/rare_diseases/data/pivus/bams/"$1}' |
	parallel --jobs 2 --col-sep "\t" "${STAR_script} {1} {2} {3} /srv/scratch/restricted/rare_diseases/data/mapping/STAR_INDEX_OVERHANG_75/ 75"

wait

STAR --genomeDir /srv/scratch/restricted/rare_diseases/data/mapping/STAR_INDEX_OVERHANG_75 --genomeLoad Remove

echo "mapping finished"
echo ""


echo "remove fastq files"
echo ""
rm $FASTQ_DIR/*merge_R1.trimmed.fastq.gz
rm $FASTQ_DIR/*merge_R2.trimmed.fastq.gz


echo "start Bam filter for mapping quality and PCR duplicates"
echo ""
cd $BAM_DIR


ls ${BAM_DIR}/*.Aligned.toTranscriptome.out.bam | parallel --jobs 10 --col-sep "\t" "${FILTER_script} {1}"

echo "Filter done"
echo ""



echo "Order by read name"
echo ""
ls *.Aligned.toTranscriptome.out_mapq30_sorted_dedup.bam | sed 's/\./\t/'|awk '{print $1}' |awk 'BEGIN{OFS="\t"}{print "/srv/scratch/restricted/rare_diseases/data/pivus/bams/"$1".Aligned.toTranscriptome.out_mapq30_sorted_dedup.bam","/srv/scratch/restricted/rare_diseases/data/pivus/bams/"$1".Aligned.toTranscriptome.out_mapq30_sorted_dedup_byread"}' |    parallel  --jobs 10 --col-sep "\t" "${SORT_SCRIPT} {1} {2}"

echo "Remove bam files"
echo ""

rm *.Aligned.toTranscriptome.out.bam
rm *.Aligned.sortedByCoord.out.bam


echo "start quantifying pivus"
cd $BAM_DIR

ls *.Aligned.toTranscriptome.out_mapq30_sorted_dedup_byread.bam | sed 's/\./\t/'|awk '{print $1}' |awk 'BEGIN{OFS="\t"}{print "/srv/scratch/restricted/rare_diseases/data/pivus/bams/"$1".Aligned.toTranscriptome.out_mapq30_sorted_dedup_byread.bam","/srv/scratch/restricted/rare_diseases/data/pivus/rsem", $1}' |    parallel   --jobs 10 --col-sep "\t" "${RSEM_script} {1} {2} {3}"



echo "Remove filtered bam files"
echo ""

rm *.Aligned.out_mapq30_sorted_dedup.ba*
rm *.Aligned.toTranscriptome.out_mapq30_sorted_dedup_byread.ba*


echo "finished quantifying pivus"
echo ""
date

