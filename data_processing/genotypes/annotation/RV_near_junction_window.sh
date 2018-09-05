#!/bin/bash
#LF
#Feb 2018

# filter for rare variants near annotated splicing junctions.

njobs=10

#bed_files_dir contains annotated bed transformed variant data
#RV_directory is directory used for output
sample_script=process_sample_RV_annotated_junctions_window.sh

date
echo "Filter each sample for RV within 20bp of for annotated junctions"
echo ""
# For each sample, filter outlier results and for rare variant nearby
ls ${bed_files_dir}/* | parallel -j $njobs --col-sep "\t" "bash ${sample_script} {1} "

wait

echo "All done!"
date
