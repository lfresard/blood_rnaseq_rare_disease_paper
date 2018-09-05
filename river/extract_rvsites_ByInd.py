#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
: Author - Yungil Kim 
: E-mail - ipw012@gmail.com

: Description

This script reads targeted regions per gene (chrom, start, stop, gencode_idx, gene_name, TSS). For each gene, search the candidate sites having MAF < 0.01 for GTEx data
and EUR population from 1k genome by navigating AF vcf files. Then, for each inidividual, write a text file having targeted sites (chrom, pos, ref, alt, gencode_idx, gene_name, major==ref?, # var, DistTSS)

: Modified by Xin, Mar 12, 2018
Xin: now output additionally: GTEx AF, 1TG AF, and gencode_idx = 0 not supplied
Xin: check REF consistency, check POS consistency, when annotating with AF.frq files
Xin: GTEx split multi-allelic sites into separate lines, tricky to handle, AF_REF is not correct
(chrom, start, ref, alt, r_allele, gene_name, major==ref?, # var, DistTSS, AF, AF_1KG, orginalGT, indID)



cd
source ./.bashrc
cat ${RAREVARDIR}/RIVER/data/rvsite/region.tss10k.txt | ${RAREVARDIR}/RIVER/extract_rvsites_ByInd.py -n 1 --id $ID --WGSvcf_in ${RAREVARDIR}/data/wgs/a_filtered_and_compressed_GTEx_WGS_vcf_file --GTExvcf_in ${RAREVARDIR}/data/wgs/a_compressed_GTEx_allele_frequency_file --EURvcf_in ${RAREVARDIR}/data/wgs/1KG/a_compressed_EUR_allele_freq_vcf_file --site_out ${RAREVARDIR}/RIVER/data/score/indiv/${ID}.txt

"""

import sys, os, re
import pysam
import numpy as np
import gzip
from optparse import OptionParser
import time
start_time = time.time()

parser = OptionParser()
parser.add_option("-n", type="int", dest="num")
parser.add_option("--id",action="store", type="string", dest="ind_id")
parser.add_option("--WGSvcf_in", dest="WGSvcf_in")    # WGS data
parser.add_option("--AFvcf_in", dest="GTExvcf_in")    # within-cohort allele frequencies
parser.add_option("--AFpop_in", dest="EURvcf_in")	# population allele frequency like 1000 genomes
parser.add_option("--site_out", dest="site_out")  # chrom, pos, ref, alt, gencode_idx, gene_name, major==ref?, # var, DistTSS
parser.add_option("--AFcutoff", type="float", dest="AFcutoff") # allele frquency cutoff within cohort samples
parser.add_option("--NONREFonly",type="string", dest="NONREFonly") # whether or not considering nonref alleles only

(options, args) = parser.parse_args()

## uploading input files
if os.path.exists(options.WGSvcf_in) and os.path.exists(options.WGSvcf_in+".tbi"):
    WGSTabix = pysam.Tabixfile(options.WGSvcf_in,'r') # WGS.vcf
if os.path.exists(options.GTExvcf_in) and os.path.exists(options.GTExvcf_in+".tbi"):
    GTExTabix = pysam.Tabixfile(options.GTExvcf_in,'r') # GTEx.vcf
if os.path.exists(options.EURvcf_in) and os.path.exists(options.EURvcf_in+".tbi"):
    EURTabix = pysam.Tabixfile(options.EURvcf_in,'r') # EUR.vcf


WGS_header = []
with gzip.open(options.WGSvcf_in,"r") as f:
    for line in f:
        eachcol = line.rstrip().split('\t')
        if eachcol[0].startswith("#CHROM") == 1:
            for colnames in eachcol:
                WGS_header.append(colnames)
       	    break



dic_WGS = {}
for indivs in WGS_header:
    dic_WGS[indivs] = WGS_header.index(indivs)

if options.AFcutoff is not None:
    rv_af_cutoff = float(options.AFcutoff)
else:
    rv_af_cutoff = 0.01

print "allele frquency cutoff: ", str(rv_af_cutoff)

# use nonref only to speed up processing
if options.NONREFonly is not None:
    NONREFonly = (options.NONREFonly in ["true", "True"])
else:
    NONREFonly = False
if NONREFonly:
    print "considering nonref alleles only"



site_out = open(options.site_out,'w')

for target_region in sys.stdin: # targeted regions per gene
    region_info = target_region.rstrip().split('\t')
    gencode_idx = 0 #region_info[3], skipped
    gene_name = region_info[3]
    tss_pos = region_info[4]
    # print region_info # timing output


    for poss_var in (WGSTabix.fetch(region_info[0],int(region_info[1]),int(region_info[2]))):


        fields_WGSvcf = poss_var.rstrip().split('\t')
	# print fields_WGSvcf
        chrom = fields_WGSvcf[0]
        pos = int(fields_WGSvcf[1])
	ref = fields_WGSvcf[3]
	alt = fields_WGSvcf[4]
        set_allele = [ref]   # ref
        set_allele.extend(re.split(',',alt)) # alt

        tempFORMAT = re.split(':', fields_WGSvcf[int(dic_WGS[options.ind_id])])
	gt_ind = tempFORMAT[0]
        gt_allele = re.split('[|/]',gt_ind) # individual-specific indices of alleles based on ref and alt

	if NONREFonly and (gt_allele[0] == "0") and (gt_allele[1] == "0"):
	    continue
	

	# a collection of all the rare alleles in the .frq file
	dic_allele_idx = {}; dic_allele_af = {}
        for target_site in GTExTabix.fetch(chrom,pos-1,pos): # AFs in GTEx European subjects
            fields_gtex = target_site.rstrip().split('\t')   
            count_allele = 0; rv_allele = [];
	    # print fields_gtex
	    if int(fields_gtex[1]) != pos: # make sure imported locus matches, tabix fetch works weirdly
		# print target_site, chrom, pos, set_allele
		continue 
	    temp_allele = fields_gtex[4].split(':')
	    if temp_allele[0] != ref:
		continue # make sure ref matches, indel and snp may appear on different lines
            for gtex_allele in fields_gtex[4:]: # summarize all rare alleles at this line of .frq
                temp_allele = gtex_allele.split(':')
		dic_allele_idx[temp_allele[0]] = count_allele; count_allele += 1
                dic_allele_af[temp_allele[0]] = temp_allele[1]
                if (float(temp_allele[1]) < rv_af_cutoff) and (float(temp_allele[1]) > 0): # GTEx rv
                    rv_allele.append(temp_allele[0])
	    # multiple lines in .frq may be summed here	
	
	# there are cases where different alt alleles are separated into different lines in vcf
	# in those cases genotypes are not correctly coded by .,0,1 only in GTEx vcf
	# also the ref frequency recorded in .frq by vcf is not correct

        if len(rv_allele) > 0:
            major_allele = dic_allele_idx.keys()[dic_allele_af.values().index(sorted(dic_allele_af.values())[-1])]

	    if major_allele not in set_allele:
		print chrom, pos, set_allele, dic_allele_idx.keys()
		print target_site

            if set_allele.index(major_allele) == 0: 
                ind_major_ref = 1; # indicator that ref. allele is major
            else:
                ind_major_ref = 0
	    
	    # a map, mark the number of occurences of rare alleles
            target_rv = {}
            for gt in gt_allele[0:2]:
		# can have two different rare alleles
                if gt == '.': continue
                elif set_allele[int(gt)] in rv_allele: # 
                    if len(target_rv) == 0:
                        target_rv[set_allele[int(gt)]] = 1
                    elif (len(target_rv) == 1) and (target_rv.keys()[0] == set_allele[int(gt)]):
                        target_rv[set_allele[int(gt)]] += 1
                    elif (len(target_rv) == 1) and (target_rv.keys()[0] != set_allele[int(gt)]):    
                        target_rv[set_allele[int(gt)]] = 1
	    
	    # count_var: how many times this rare variant appear in this genotype
	    # different rv_allele will be output on different lines
	    # same rv_allele will be single line with count_var indicating occurrences 
            for ind_rv,count_var in target_rv.iteritems():
                fields_gen1k_var = []; out_site = []; find_frq = False
                for gen1k_var in EURTabix.fetch(chrom,pos-1,pos):
                    fields_gen1k_var = gen1k_var.rstrip().split('\t')
		    if int(fields_gen1k_var[1]) != pos:
			continue # make sure the pos matches the 1KG.frq entry

		    # if only the ref matches
		    temp_allele = fields_gen1k_var[4].split(':')
		    if temp_allele[0] != ref:
			continue
		    else:
			find_frq = True

		    for gen1k_allele in fields_gen1k_var[4:]:
			temp_allele = gen1k_allele.split(':')
              		if (temp_allele[0] == ind_rv) and (float(temp_allele[1]) < 0.01):
                            out_site = [str(chrom),str(pos),str(ref),str(alt),str(ind_rv),str(gene_name), \
                                    	str(ind_major_ref),str(count_var),str(abs(int(tss_pos)-int(pos))), \
					str(dic_allele_af[ind_rv]), str(temp_allele[1]), \
					gt_ind, options.ind_id]
			    break

                
                if not find_frq: # no variant in EUR population
                    out_site = [str(chrom),str(pos),str(ref),str(alt),str(ind_rv),str(gene_name), \
                            str(ind_major_ref),str(count_var),str(abs(int(tss_pos)-int(pos))), \
			    str(dic_allele_af[ind_rv]), str(0.0), \
			    gt_ind, options.ind_id]
                
                if len(out_site) > 0:
                    site_out.write("\t".join(out_site)+"\n")
		    # print out_site
		    

        elif len(rv_allele) == 0:
            continue

#print("--- %s hours ---" % ((time.time() - start_time)/3600))

#f = open('${RAREVARDIR}/RIVER/data/score/indiv/done.txt', 'a+')
#f.write("%s\n" % (str(options.num)))
#f.close()






