# Expression outlier scripts

Scripts to re-create the expression outlier enrichment analyses as displayed in plots 2A,B and S5C,D,E and expression outlier counts as displayed in plot 2C and S7

## Expression outlier enrichment

### Steps

##### File: `logistic_percentile.r`

* Define threshold for global expression outliers (for the paper, the threshold was outlier in >100 genes at abs(Z-score)>=4)
* Read in corrected gene expression count data (`corrected_counts/gene/blood/[gene_count].txt`) and sample metadata (`[metadata_file].txt`)
* Subset samples to include blood samples only
* Center and scale gene expression counts matrix to generate Z-scores
* Read in ExAC gene constraint metrics `[ExAC_constraint_metrics].txt` and add ensgene ID
* Find global expression outlier samples and remove from gene expression counts matrix
* Define percentile thresholds for expression outlier calling
* For each percentile threshold, count case/control expression outliers and ExAC constraint metrics for each gene and run logistic model for each mutation class
* Write analysis to file
* Write global expression outliers to file

### Additional Steps 

Scripts that repeat the steps carried out by `logistic_percentile.r`, but with added functionality

##### File: `logistic_percentile_vary_control_n.r`

This script is intended to be run in parallel, using a different random subsets of external control samples in each run. The script was designed to be submitted using the SLURM job management system using the command `Rscript scg_super_sva_logistic_vary_global.R ${SLURM_ARRAY_TASK_ID}` where `${SLURM_ARRAY_TASK_ID}` is an integer value between 1:100000 for each of the 100000 permutations performed. 

* Repeat pipeline for different external control sample thresholds defined in `n_dgn_control`
* Perform SVA correction using uncorrected gene count data (`[raw_gene_counts].txt`)
* Read a global expression outlier vector (`[global_out[.txt`) output by `logistic_percentile.r`
* Read a file containing 10000 columns of randomly chosen integers from 1 to the total number of external control samples (the column chosen on each run is dependent on current integer value of `${SLURM_ARRAY_TASK_ID}` for a given run)
* Repeat logistic model pipeline for rare disease case and unaffected family member control label switching
* Write output

##### File: `logistic_permutecasecontrol_parallel.r`

This script...

## Expression outlier counts

### Steps
* 