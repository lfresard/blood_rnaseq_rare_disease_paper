# Surrogate Variable Analysis (SVA) data correction pipeline

Script to correct expression data using surrogate variables + regression splines

### Steps

##### File: `sva_pipeline.r`

* Read in raw count data output by RSEM
* Filter genes for minimum expression threshold (in the paper - TPM>0.5 in >50% of samples from each batch)
* Log transform
* Scale and center to generate expression Z-scores
* Run SVA
* Regress out SVs and SV+regression splines for SVs significantly associated with sample batch
* Run QC checks