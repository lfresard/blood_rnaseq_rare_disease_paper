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
bash splicing_outlier_analysis.sh \
	metadatafile.tsv \ # metadata file to use
	Blood  \ # tissue of analysis
	$junc_dir \ # folder containing filtered junction files
	$output_dir \ # output folder
	TRUE \ # wether or not to perform on samples from the freeze only
	FALSE \ # wether or not to include DGN samples
	TRUE > log_file.txt 2>&1 & # wether or not to include PIVUS in the analysis
```
## 2- Generating Z-scores from the ratios
This script impute missing splicing ratios and get Z-scores for splicing data.
```
Rscript splicing_ratio_to_zscores.R
```

## 3- Cross Z-score results with RV information
Need to sort the outlier file first

* distance=20

```
sort -k1,1 -k2,2n  splicing_outlier_file.txt>  splicing_outlier_file_sorted.txt

bash splicing_outlier_genes_with_RV_window.sh  splicing_outlier_file_sorted.txt > <log_file> 2>&1 &
```
