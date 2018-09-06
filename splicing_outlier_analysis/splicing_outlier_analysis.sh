#!/bin/bash

# This script handle the splice junction analysis for every rare disease case on a given tissue
# start point is metadata file
# gives the list of case and controls for each case.

#For each case we need to:
#-create a folder
#-make a list of case and controls
#-run the splicing analysis script
#-filter results.


##--- Variables
metadata_file=$1 # file containing data origin and status information 
tissue=$2 # what tissue is the analysis made (example "Blood")
junc_dir=$3 # directory containg filtered junctions
out_dir=$4 # output directory
freeze=$5 # set of samples to use (yes or no)
DGN=$6 # use DGN in control (yes or no)
PIVUS=$7 # use PIVUS as controls (yes or no)


date
echo ""
# Create output dir if does not exist
if [[ ! -e $out_dir ]]; then
    mkdir $out_dir
elif [[ ! -d $out_dir ]]; then
    echo "$out_dir already exists but is not a directory" 1>&2
fi


# Organize results directory and create case control lists
echo "Organize results directory and create case control lists"
echo ""
/users/lfresard/R_3.3.2/bin/Rscript make_list_junction_files.R --meta $metadata_file --tissue $tissue --outdir $out_dir --juncdir $junc_dir --freeze $freeze --DGN $DGN --PIVUS $PIVUS


wait

# Run junction outlier analysis
echo "Run junction outlier analysis"
echo ""

if [[ $DGN == TRUE ]]; then
	outPrefix='RD_DGN_freeze'
	list_junc=list_junctions_freeze_RD_DGN.txt
	meta_file=sample_affected_status_freeze_RD_DGN.tsv

elif [[ $PIVUS == TRUE ]]; then 
	outPrefix='RD_PIVUS_freeze'
	list_junc=list_junctions_freeze_RD_PIVUS.txt
	meta_file=sample_affected_status_freeze_RD_PIVUS.tsv

else
	outPrefix='RD_freeze'
	list_junc=list_junctions_freeze_RD.txt
	meta_file=sample_affected_status_freeze_RD.tsv
fi 

python splicing_outlier.py --juncfiles ${list_junc} --outprefix ${outPrefix} --meta_file ${meta_file} --rundir ${out_dir}
wait


echo "Done"
echo ""

date
