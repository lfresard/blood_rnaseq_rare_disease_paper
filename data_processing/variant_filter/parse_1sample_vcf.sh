#!/bin/bash

vcf_file=$1

# processes vcf with one sample

# accomodate following FORMAT
# CGS	GT:AD:DP:GQ:PL
# CGS	GT:GQ:PL
# CHEO	GT
# CHEO	GT:AD:DP:GQ:PL
# CHEO	GT:GQ:PL
# CHEO	GT:PL
# CHEO	GT:PL:DP:SP:GQ
# UDN	GT:AD
# UDN	GT:AD:DP
# UDN	GT:AD:DP:GQ:PGT:PID:PL
# UDN	GT:AD:DP:GQ:PL
# UDN	GT:AD:DP:PGT:PID
# UDN	GT:AD:GQ:PGT:PID:PL
# UDN	GT:AD:GQ:PL
# UDN	GT:AD:PGT:PID
# UDN	GT:VR:RR:DP:GQ

# AD field cannot be normed by bcftools norm --multiallelics -any
# need to process manually

# make sure all fields filled with "NA" if not available

printf "sample_id\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tGT\tDP\tGQ\tPL[0]\tPL_GT\tAD[0]\tAD_H1\tAD_H2\tH1\tH2\tvflag\n"

if [[ ${vcf_file} =~ (RD[0-9]+).vcf.gz ]]; then
	id=${BASH_REMATCH[1]};
	zcat ${vcf_file} | \
		awk -F'\t' -vID=$id \
		'BEGIN{OFS="\t"} \
                $2 ~ /[0-9]+/ { \
			FT=$7; \
                        PL[1]="NA";PL[2]="NA";\
                        AD[1]="NA";AD[2]="NA";AD[3]="NA";\
                        GT="NA";GQ="NA";DP="NA"; \
			H1=0;H2=0; \
			H1a="";H2a=""; \
			tempGT[1]="NA";tempGT[2]="NA"; \
			tempAD[1]="NA";tempAD[2]="NA";tempAD[3]="NA"; \
                        split($9, FORMAT, ":"); \
                        split($10, FIELD, ":"); \
			split($4","$5, ALT, ","); \
                        for(i=1;i<=length(FORMAT);i++){ \
                                if(FORMAT[i]=="DP"){DP=FIELD[i];}; \
                                if(FORMAT[i]=="GT"){GT=FIELD[i]; split(GT,tempGT,/[/|]/); if(tempGT[1] ~ /[0-9]+/){H1=tempGT[1];H1a=ALT[H1+1]}; if(tempGT[2] ~ /[0-9]+/){H2=tempGT[2];H2a=ALT[H2+1]};}; \
                                if(FORMAT[i]=="GQ"){GQ=FIELD[i];}; \
                                if(FORMAT[i]=="PL"){split(FIELD[i], tempPL, ","); if(H1>=H2){o=(H1+1)*H1/2+H2}else{o=(H2+1)*H2/2+H1;}; PL[1]=tempPL[1]; PL[2]=tempPL[o+1]}; \
                                if(FORMAT[i]=="AD"){split(FIELD[i], tempAD, ","); AD[1]=tempAD[1]; AD[2]=tempAD[H1+1]; AD[3]=tempAD[H2+1]}; \
				if(FORMAT[i]=="RR"){tempAD[1]=FIELD[i]; AD[1]=tempAD[1]; AD[2]=tempAD[H1+1]; AD[3]=tempAD[H2+1]}; \
				if(FORMAT[i]=="VR"){tempAD[2]=FIELD[i]; AD[1]=tempAD[1]; AD[2]=tempAD[H1+1]; AD[3]=tempAD[H2+1]}; \
                        }; \
			bitF = (FT == "NA" || FT=="." || FT == "num prev seen samples > 30" || FT=="PASS"); \
			bitGQ = (GQ == "NA" || GQ == "." || GQ >= 20); \
			bitDP = (DP == "NA" || DP == "." || DP >= 20); \
			bitPL = (PL[1] == "NA" || PL[2] < 20); \
			if( AD[1] == "NA" || DP == "NA" || DP == 0 ){bitAD = 1}else{ \
				if( H1 == H2 ){bitAD = (AD[2]/DP > 0.8)}else{ \
					bitAD = (AD[2]/DP > 0.2 && AD[3]/DP > 0.2 && (AD[2]+AD[3])/DP > 0.8) }}\
			if( bitF && bitGQ && bitDP && bitPL && bitAD ){flag = "PASS";}else{flag = "FAIL";}\
			print ID,$1,$2,$3,$4,$5,$6,FT, \
                	GT,DP,GQ,PL[1],PL[2],AD[1],AD[2],AD[3],H1":"H1a,H2":"H2a,flag;}'
fi

# need to incorporate variant quality filter extraction into river pipeline
# in case of VR:RR notation, multi-allelic sites are ambiguous

# awk variables are dynamically typed, numeric looking string will be converted automatically

