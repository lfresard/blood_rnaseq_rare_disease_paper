# Identify ASE outliers and assess association with phenotype via Amelie

From ASE output to z-scores per gene and Amelie scores

### Steps
* Collapse ASE data across all GTEx v7 tissues per individual per gene using same filters as RDS (read count > 20 and remove absolutes) and selecting the max abs(0.5 - refRatio) value across tissues per individual per gene
* Merge GTEx and RDS ASE data and filter to sites and genes seen in both
* For each gene, scale the reference ratios across all sites within the gene across RDS and GTEx
* Call an individual-gene an ASE outlier if it contains the most extreme up or down site for that gene across all samples
* Remove global outlier sample
* Assess RDS case outlier genes for association with that case's HPO terms via Amelie
* Over 100 random iterations, select a new random set of genes per case and calculate phenotype association via Amelie
* Compare distribution of matched and random Amelie scores

### Scripts
collapse_GTEx_ASE.R - requires all GTEx v7 ASE files per tissue and outputs, RDS ase data and sample metadata and outputs filtered and collapsed GTEx data in `gtexV7_ase_filteredSitesForRDS_maxPerTissue.txt` for use downstream

get_ASE_outliers_withGTEx.R - requires `gtexV7_ase_filteredSitesForRDS_maxPerTissue.txt` from above and RDS ase data, merges the data frames and scales reference ratios to obtain z-scores, outputting those values to `RDS_GTEX_ASE_zscores.txt` and outlier counts to `RDS_ASE_outlier_counts.txt`

permute_amelie_calls_ASEoutliers.R - Reads in z-scores from above, as well as RDS metadata and HPO information and runs Amelie to obtain phenotype association scores for all ASE outlier genes, and then runs 100 iterations of selecting random genes to score and produces a plot of the distribution of scores when genes are matched vs random

```
Rscript collapse_GTEx_ASE.R --asedir [ase directory] --datadir [data output directory] --rdsdir [RDS data directory] &
Rscript get_ASE_outliers_withGTEx.R --datadir [data output directory] &
Rscript permute_amelie_calls_ASEoutliers.R --datadir [data output directory] --ameliedir [Directory with amelie scripts] &
```
