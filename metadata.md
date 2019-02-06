# Metadata file format

This file contains the column names from the metadata file that are needed in the pipeline
* `sample_id`: identifier of the RNAseq smaple
* `affected_status`: Case/Control
* `in_freeze`: yes/no, wether sample is in teh current freeze of analysis
* `is_RD`: yes/no, wether the sample is a rare disease sample
* `source_of_RNA`: tissue of origin for the sample (i.e Blood, Fibroblast,...)
* `HPO_terms_ID`: Comma separated list of HPO ID linked to the phenotype of the sample
* `status`: PASSED/FAILED, wether the RNAseq was successful
* `variant_data`: genome/exome/none

