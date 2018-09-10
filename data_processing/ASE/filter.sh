#!/bin/bash



LOGFILE=$RARE_DIS_DIR/analysis/ase_analysis/Nikki_Nov2017/runlogs/Nov30_1

ASERC_out=$RARE_DIS_DIR/data/ase/ASEReadCounter_output
filtered_out=$RARE_DIS_DIR/data/ase/filtered

for file in $ASERC_out/RD001*.ase.csv
do
    fileh=${file##*/}
    SAMPID=$(echo $fileh | cut -f 1 -d '.')
    if [ -s $FILTERED_OUT/$SAMPID.ase.filtered.csv ]
    then
	echo "$SAMPID ase filtered already exists, skipping" `date "+%Y %m %d - %H:%M:%S"` >> $LOGFILE
    else
	echo "Filtering $SAMPID" `date "+%Y %m %d - %H:%M:%S"` >> $LOGFILE
	awk 'NR == 1 {print $0"\trefRatio"} NR > 1 && $6 > 0 && $7 > 0 {print $0"\t"$6/$8}' $file > $FILTERED_OUT/$SAMPID.ase.filtered.csv
    fi
done
