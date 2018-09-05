#LF
#Jan 2018

# This script runs RD_vcf_toRVbed.sh
# filters for Rare variants
# transforms into bed file
# get genes name overlapping with rare variant

vcf_dir=$1
cd ${vcf_dir}
nproc=10

vcf_suffix="_homogenized_gnomad_cadd.vcf.gz"

ls RD*${vcf_suffix} | parallel --jobs $nproc "bash /users/lfresard/repos/rare_disease/scripts/vcf_handling/RD_vcf_to_RVbed.sh {}"
