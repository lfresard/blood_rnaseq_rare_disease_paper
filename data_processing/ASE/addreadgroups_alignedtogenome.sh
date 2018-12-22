#!/bin/bash

OUTBAM=$RARE_DIS_DIR/data/filtered_bam/for_freeze/AlignedToGenomeReadGroups

for file in $RARE_DIS_DIR/data/filtered_bam/for_freeze/AlignedToGenome/*.Aligned.sortedByCoord.out_mapq30_sorted_dedup.bam
do
    fileh=${file##*/}
    filen=${fileh%.*}
    echo $filen
    java -jar /software/picard-tools/1.92/AddOrReplaceReadGroups.jar \
	 I=$file \
	 O=$OUTBAM/$filen.bam \
	 RGID=1 \
	 RGLB=lib1 \
	 RGPL=illumina \
	 RGPU=unit1 \
	 RGSM=20

    #also index the file
    /software/samtools/samtools-1.2/bin/samtools index $OUTBAM/$filen.bam
done
