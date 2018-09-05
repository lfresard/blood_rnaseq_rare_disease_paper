#!/bin/bash

RUNNAME=$1

LOGFILE=/srv/scratch/restricted/rare_diseases/analysis/ase_analysis/Nikki_Nov2017/runlogs/$RUNNAME

BAMPATH=/srv/scratch/restricted/rare_diseases/data/filtered_bam/for_freeze/AlignedToGenomeReadGroups
WORKINGDIR=/srv/scratch/restricted/rare_diseases/data/ase/$RUNNAME
GATKPATH=/srv/scratch/restricted/rare_diseases/analysis/ase_analysis/GATK
GENOME=/mnt/lab_data/montgomery/shared/genomes/hg19/hg19.fa
VCFFILEPATH=/srv/scratch/restricted/rare_diseases/data/vcfs/intermediates/for_freeze
SCRIPTS=/srv/scratch/restricted/rare_diseases/analysis/ase_analysis/Nikki_Nov2017/scripts

mkdir $WORKINGDIR
mkdir $WORKINGDIR/VCF
mkdir $WORKINGDIR/ASERC_out
mkdir $WORKINGDIR/filtered
mkdir $WORKINGDIR/stats
mkdir $WORKINGDIR/gene_table

for VCFFILE in $VCFFILEPATH/RD*.vcf.gz
do
    fileh=${VCFFILE##*/}
    SAMPID=$(echo $fileh | cut -f 1 -d '_') ##files are in format ID_interim.vcf.gz

    #if I haven't processed it yet
    if [ -s $OUTPATH/$SAMPID.ase.csv ]
    then
	echo "$SAMPID ase already exists, skipping" `date "+%Y %m %d - %H:%M:%S"` >> $LOGFILE
    else

	echo -n "$SAMPID started:"`date "+%Y/%m/%d - %H:%M:%S"` >> $LOGFILE

	#merge biallelic sites into multiallelic
	#unzip the file, else aserc gives an empty output
	gzip -dc < $VCFFILE | sed 's/ID=GT,Number=1,Type=Integer,/ID=GT,Number=1,Type=String,/g' | bcftools call --skip-variants indels - > $WORKINGDIR/VCF/$SAMPID.vcf

	gzip -dc < $VCFFILE | sed 's/ID=GT,Number=1,Type=Integer,/ID=GT,Number=1,Type=String,/g' | vcftools --vcf - --remove-indels --recode --recode-INFO-all --stdout | bcftools norm --multiallelics +snps --no-version -Ov - > $WORKINGDIR/VCF/$SAMPID.filtered.vcf

	awk '{ if ($1 ~ /[#]/) print $0; else exit;}' $WORKINGDIR/VCF/$SAMPID.filtered.vcf > $WORKINGDIR/VCF/$SAMPID.filtered.hg19sort.vcf

	awk -F $'\t' 'BEGIN {OFS = FS} { if ($1 ~ /[#]/); else if ($1 == "19") print 20.6,$0; else if ($1 == "22") print 20.7,$0; else if ($1 ~/[1-9]/) print $1,$0; else if ($1 == "X") print 7.5,$0; else if ($1 == "Y") print 20.5,$0;  }' $WORKINGDIR/VCF/$SAMPID.filtered.vcf | sort -k1,1n -k3,3n | awk '{for (i=2; i<NF; i++) printf $i "\t"; print $NF}' | awk '{ if($0 !~ /^#/) print "chr"$0; else if(match($0,/(##contig=<ID=)(.*)/,m)) print m[1]"chr"m[2]; else print $0 }' >> $WORKINGDIR/VCF/$SAMPID.filtered.hg19sort.vcf

	echo -n "\t Running ASE"`date "+%Y/%m/%d - %H:%M:%S"` >> $LOGFILE

	java -jar $GATKPATH/GenomeAnalysisTK.jar \
	     -R $GENOME \
	     -T ASEReadCounter \
	     -o $WORKINGDIR/ASERC_out/$SAMPID.ase.tsv \
	     -I $BAMPATH/$SAMPID.Aligned.sortedByCoord.out_mapq30_sorted_dedup.bam \
	     -sites $WORKINGDIR/VCF/$SAMPID.filtered.hg19sort.vcf \
	     -U ALLOW_N_CIGAR_READS \
	     -minDepth 10 \
	     --minMappingQuality 10 \
	     --minBaseQuality 2 \
	     -drf DuplicateRead

	#filter reads
	echo -n "\t Post-processing"`date "+%Y/%m/%d - %H:%M:%S"` >> $LOGFILE
	awk 'NR == 1 {print $0"\trefRatio"} NR > 1 && $6 > 0 && $7 > 0 {print $0"\t"$6/$8}' $WORKINGDIR/ASERC_out/$SAMPID.ase.tsv > $WORKINGDIR/filtered/$SAMPID.ase.filtered.tsv

	python $SCRIPTS/ase_stats.py $WORKINGDIR/filtered/$SAMPID.ase.filtered.tsv $WORKINGDIR/stats/$SAMPID.ase.stats.tsv

	# contigpositionvariantIDrefAllelealtAllelerefCountaltCounttotalCountlowMAPQDepthlowBaseQDepth11.rawDepthotherBasesimproperPairsrefRatio15.eb_ref_ratioalpha1beta1eb_loweb_highPEP21qvalue >>>> contig pos-1 pos altAllele,rawDepth,eb_ref_ratio,qvalue
	  awk 'NR > 1 {print $1"\t"$2-1"\t"$2"\t"$7" "$11" "$15" "$21}' $WORKINGDIR/stats/$SAMPID.ase.stats.tsv | bedtools intersect -wa -wb -a stdin -b $SCRIPTS/UCSC_geneexon.bed | awk 'BEGIN{OFS="\t"; print "#Chr\tPOS\tGene\tExonCount\tAltReads\tTotalReads\tebRefRatio\tQValue"} {print $1,$3,$11,$12,$4,$5,$6,$7}' > $WORKINGDIR/gene_table/$SAMPID.ase.genetable.tsv

	echo "\t Done!"`date "+%Y/%m/%d - %H:%M:%S"` >> $LOGFILE
    fi
done


echo "#!/bin/bash

RUNNAME=$1

LOGFILE=/srv/scratch/restricted/rare_diseases/analysis/ase_analysis/Nikki_Nov2017/runlogs/$RUNNAME

BAMPATH=/srv/scratch/restricted/rare_diseases/data/filtered_bam/for_freeze/AlignedToGenomeReadGroups
WORKINGDIR=/srv/scratch/restricted/rare_diseases/data/ase/$RUNNAME
GATKPATH=/srv/scratch/restricted/rare_diseases/analysis/ase_analysis/GATK
GENOME=/mnt/lab_data/montgomery/shared/genomes/hg19/hg19.fa
VCFFILEPATH=/srv/scratch/restricted/rare_diseases/data/vcfs/intermediates/for_freeze
SCRIPTS=/srv/scratch/restricted/rare_diseases/analysis/ase_analysis/Nikki_Nov2017/scripts

mkdir $WORKINGDIR
mkdir $WORKINGDIR/VCF
mkdir $WORKINGDIR/ASERC_out
mkdir $WORKINGDIR/filtered
mkdir $WORKINGDIR/stats
mkdir $WORKINGDIR/gene_table

for VCFFILE in $VCFFILEPATH/RD*.vcf.gz
do
    fileh=${VCFFILE##*/}
    SAMPID=$(echo $fileh | cut -f 1 -d '_') ##files are in format ID_interim.vcf.gz

    #if I haven't processed it yet
    if [ -s $OUTPATH/$SAMPID.ase.csv ]
    then
        echo "$SAMPID ase already exists, skipping" `date "+%Y %m %d - %H:%M:%S"` >> $LOGFILE
    else

echo -n "$SAMPID started:"`date "+%Y/%m/%d - %H:%M:%S"` >> $LOGFILE

#merge biallelic sites into multiallelic
#unzip the file, else aserc gives an empty output
gzip -dc < $VCFFILE | sed 's/ID=GT,Number=1,Type=Integer,/ID=GT,Number=1,Type=String,/g' | bcftools call --skip-variants indels - > $WORKINGDIR/VCF/$SAMPID.vcf

gzip -dc < $VCFFILE | sed 's/ID=GT,Number=1,Type=Integer,/ID=GT,Number=1,Type=String,/g' | vcftools --vcf - --remove-indels --recode --recode-INFO-all --stdout | bcftools norm --multiallelics +snps --no-version -Ov - > $WORKINGDIR/VCF/$SAMPID.filtered.vcf 

awk '{ if ($1 ~ /[#]/) print $0; else exit;}' $WORKINGDIR/VCF/$SAMPID.filtered.vcf > $WORKINGDIR/VCF/$SAMPID.filtered.hg19sort.vcf

awk -F $'\t' 'BEGIN {OFS = FS} { if ($1 ~ /[#]/); else if ($1 == "19") print 20.6,$0; else if ($1 == "22") print 20.7,$0; else if ($1 ~/[1-9]/) print $1,$0; else if ($1 == "X") print 7.5,$0; else if ($1 == "Y") print 20.5,$0;  }' $WORKINGDIR/VCF/$SAMPID.filtered.vcf | sort -k1,1n -k3,3n | awk '{for (i=2; i<NF; i++) printf $i "\t"; print $NF}' | awk '{ if($0 !~ /^#/) print "chr"$0; else if(match($0,/(##contig=<ID=)(.*)/,m)) print m[1]"chr"m[2]; else print $0 }' >> $WORKINGDIR/VCF/$SAMPID.filtered.hg19sort.vcf

echo -n "\t Running ASE"`date "+%Y/%m/%d - %H:%M:%S"` >> $LOGFILE

java -jar $GATKPATH/GenomeAnalysisTK.jar \
     -R $GENOME \
     -T ASEReadCounter \
     -o $WORKINGDIR/ASERC_out/$SAMPID.ase.tsv \
     -I $BAMPATH/$SAMPID.Aligned.sortedByCoord.out_mapq30_sorted_dedup.bam \
     -sites $WORKINGDIR/VCF/$SAMPID.filtered.hg19sort.vcf \
     -U ALLOW_N_CIGAR_READS \
     -minDepth 10 \
     --minMappingQuality 10 \
     --minBaseQuality 2 \
     -drf DuplicateRead

#filter reads
echo -n "\t Post-processing"`date "+%Y/%m/%d - %H:%M:%S"` >> $LOGFILE
awk 'NR == 1 {print $0"\trefRatio"} NR > 1 && $6 > 0 && $7 > 0 {print $0"\t"$6/$8}' $WORKINGDIR/ASERC_out/$SAMPID.ase.tsv > $WORKINGDIR/filtered/$SAMPID.ase.filtered.tsv

python $SCRIPTS/ase_stats.py $WORKINGDIR/filtered/$SAMPID.ase.filtered.tsv $WORKINGDIR/stats/$SAMPID.ase.stats.tsv

# contigpositionvariantIDrefAllelealtAllelerefCountaltCounttotalCountlowMAPQDepthlowBaseQDepth11.rawDepthotherBasesimproperPairsrefRatio15.eb_ref_ratioalpha1beta1eb_loweb_highPEP21qvalue >>>> contig pos-1 pos altAllele,rawDepth,eb_ref_ratio,qvalue
  awk 'NR > 1 {print $1"\t"$2-1"\t"$2"\t"$7" "$11" "$15" "$21}' $WORKINGDIR/stats/$SAMPID.ase.stats.tsv | bedtools intersect -wa -wb -a stdin -b $SCRIPTS/UCSC_geneexon.bed | awk 'BEGIN{OFS="\t"; print "#Chr\tPOS\tGene\tExonCount\tAltReads\tTotalReads\tebRefRatio\tQValue"} {print $1,$3,$11,$12,$4,$5,$6,$7}' > $WORKINGDIR/gene_table/$SAMPID.ase.genetable.tsv

echo "\t Done!"`date "+%Y/%m/%d - %H:%M:%S"` >> $LOGFILE
fi
done" > $LOGFILE


