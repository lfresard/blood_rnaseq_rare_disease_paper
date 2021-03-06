#!/bin/bash
#LF
#STAR_2.4.0j
set -o nounset -o pipefail

GENOME_FASTA=hg19.fa
GTF_FILE=gencode.v19.annotation.gtf

# read command line arguments: 
read1=$1
read2=$2
prefix=$3
INDEX_DIRECTORY=$4 
OVERHANG=$5

cmd=\
"STAR \
    --genomeDir $INDEX_DIRECTORY \
    --genomeFastaFiles $GENOME_FASTA \
    --sjdbGTFfile $GTF_FILE \
    --sjdbOverhang $OVERHANG \
    --readFilesIn $read1 $read2 \
	--outFileNamePrefix ${prefix}. \
	--readFilesCommand zcat \
	--outSAMattributes NH HI AS NM MD nM\
	--outFilterType BySJout \
	--runThreadN 4 \
	--outSAMtype BAM SortedByCoordinate \
	--quantMode TranscriptomeSAM\
	--genomeLoad NoSharedMemory \
	--limitBAMsortRAM 15000000000"
    
date
echo $cmd
$cmd
date

