#!/bin/bash

## wrapper script to call amelie from R with input file and prefix as command line arguments (called from permute_amelie_calls_ASEoutliers.R)
inputfile=$1
itpre=$2
while read id genes hpo
do
	echo $id
	python amelie_per_sample.py --genes $genes --hpo $hpo --outprefix $itpre$id
done < $inputfile
