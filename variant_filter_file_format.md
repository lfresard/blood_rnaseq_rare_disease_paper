This is the format of the file created by [generate_variantfilters.sh](./data_processing/genotypes/variant_filter/generate_variantfilters.sh)


* `sample_id`: sample ID
* `CHROM`: Chromosome
* `POS`: Position
* `ID`: Identifier
* `REF`: Reference Allele
* `ALT`: Alternate Allele
* `QUAL`: Phred-scaled quality score
* `FILTER`: PASS if this position has passed all filters, i.e. a call is made at this position.
* `GT`: genotype
* `DP`: read depth at this position for this sample
* `GQ`: genotype quality, encoded as a phred quality 
* `PL[0]`: "Normalized" Phred-scaled likelihoods of being 0/0
* `PL_GT`: "Normalized" Phred-scaled likelihoods of called genotype
* `AD[0]`: unfiltered allele depth
* `AD_H1`: 
* `AD_H2`
* `H1`
* `H2`
* `vflag`: PASS or FAIL filters
