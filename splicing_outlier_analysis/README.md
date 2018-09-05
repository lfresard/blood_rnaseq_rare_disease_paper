# Perform splicing outlier analysis

This set of scripts is generating junction ratios and caluclating splicing Z-scores for every sample

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites
To get started you will need
* Metadata file
* Tissue to perform the analysis on
* Junction files generated during STAR alignment (here they have been filtered for junctions with at least 10 reads uniquely spanning)


## 1- Generating junction ratios
```
bash /users/lfresard/repos/rare_disease/scripts/splicing_analysis/outlier/splicing_outlier_analysis.sh \
	/srv/scratch/restricted/rare_diseases/data/metadata/2018_03_09_Rare_Disease_Metadata.tsv \ # metadata file to use
	Blood  \ # tissue of analysis
	/srv/scratch/restricted/rare_diseases/data/splicing/juncfiles/all_filteredjunc \ # folder containing filtered junction files
	/srv/scratch/restricted/rare_diseases/analysis/outlier_analysis/splicing/for_freeze \ # output folder
	TRUE \ # wether or not to perform on samples from the freeze only
	FALSE \ # wether or not to include DGN samples
	TRUE > log_splicing_outlier_RD_PIVUS_2018_03_13.txt 2>&1 & # wether or not to include PIVUS in the analysis
```
## 2- Generating Z-scores from the ratios
This script impute missing splicing ratios and get Z-scores for splicing data.
```
Rscript /users/lfresard/repos/rare_disease/scripts/splicing_analysis/outlier/splicing_ratio_to_zscores.R
```

## 3- Cross Z-score results with RV information
Need to sort the outlier file first

* distance=20
* outdir=/srv/scratch/restricted/rare_diseases/analysis/outlier_analysis/splicing/for_freeze/RV
* file="/srv/scratch/restricted/rare_diseases/data/bed/for_freeze/${sample}_homogenized_gnomad_cadd_RV_withgene.bed.gz"
* outlier_file="/srv/scratch/restricted/rare_diseases/analysis/outlier_analysis/splicing/for_freeze/RD_freeze_junc_outliers_PCAratio_blood_nf_PIVUS_withPIVsamples_withmisdata_sorted.txt" 
```
sort -k1,1 -k2,2n  RD_freeze_junc_outliers_PCAratio_blood_nf_PIVUS_withPIVsamples_withmisdata.txt>  RD_freeze_junc_outliers_PCAratio_blood_nf_PIVUS_withPIVsamples_withmisdata_sorted.txt

bash /users/lfresard/repos/rare_disease/scripts/splicing_analysis/outlier/splicing_outlier_genes_with_RV_window.sh  /srv/scratch/restricted/rare_diseases/analysis/outlier_analysis/splicing/for_freeze/RD_freeze_junc_outliers_PCAratio_blood_nf_PIVUS_withPIVsamples_withmisdata_sorted.txt > <log_file> 2>&1 &
```
