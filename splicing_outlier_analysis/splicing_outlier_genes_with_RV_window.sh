#!/bin/bash

# filter splicing results for genes with RV nearby
njobs=10

samples=sample_affected_status_freeze_RD.tsv
outlier_results=$1
sample_script=process_sample_RV_window.sh
outlier_RV_suffix=_outliers_RV_window.txt
outdir=<outdir>
outlier_RV_res=$(basename "$outlier_results" .txt)_RV_window.txt
cd ${outdir}

date
echo "Filter each sample for its RV"
echo ""
# For each sample, filter outlier results and for rare variant nearby
cat $samples | awk '{print $1}' |awk -v outlier_res=${outlier_results} 'BEGIN{OFS="\t"}{print $1 ,outlier_res}' | parallel -j $njobs --col-sep "\t" "bash ${sample_script} {1} {2}"

wait

echo "Concatenate results per sample"
echo ""
# Cancatenate results per sample
cat *${outlier_RV_suffix} > $outlier_RV_res
wait

# Remove intermediate files
echo "Remove intermediate files"
echo ""
rm *${outlier_RV_suffix}

echo "All done!"
date
