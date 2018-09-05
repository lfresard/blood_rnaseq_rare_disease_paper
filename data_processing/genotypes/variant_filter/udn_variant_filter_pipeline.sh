#!/bin/bash


bash generate_variantfilters.sh | parallel --jobs 20

wait

# generate_sitefilters.sh

wait

bash combine_filters.sh

wait

bash filter_allvcf.sh | parallel --jobs 20

wait

# bash attach_variantfilter.sh



