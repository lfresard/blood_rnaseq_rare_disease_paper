#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author - Yungil Kim
:E-mail - ipw012@gmail.com

:Description:

From each of individual.txt files, This script extract all of features including chromHMM, Phylop, DANN scores and combine all of features simultaneously

:Modified by Xin, Mar 12, 2018


cat ${RAREVARDIR}/RIVER/data/score/indiv/${ID}.txt |
${RAREVARDIR}/RIVER/code/extract_scores_combined.py -n $count_ind --id $ID
--af_in ${RAREVARDIR}/data/wgs/GTEx_af.vcf.gz --wgs_in
filtered_and_compressed_GTEx_WGS_vcf_file --anno_in
${RAREVARDIR}/data/wgs/GTEx_vep.vcf.gz --cadd_in
${RAREVARDIR}/RIVER/data/CADD_whole_genome_SNVs_inclAnno.tsv.gz --dann_in
${RAREVARDIR}/RIVER/data/DANN_whole_genome_SNVs.tsv.bgz --chromHMM_in
${RAREVARDIR}/RIVER/data/wgEncodeBroadHmmGm12878HMM.sorted.hg19.bed.txt.gz
--phylop_in ${RAREVARDIR}/RIVER/data/phyloP100way.txt.gz --score_out
${RAREVARDIR}/RIVER/data/score/indiv/${ID}.${count_ind}.score.nuc.txt

"""

import sys, os, re
import pysam
import numpy as np
import gzip
from optparse import OptionParser
import time

# start_time = time.time()

parser = OptionParser()

parser.add_option("-n", type="int", dest="num")
parser.add_option("--id",action="store", type="string", dest="ind_id")

parser.add_option("--af_in", dest="af_in")          # af vcf GTEX
parser.add_option("--wgs_in", dest="wgs_in")        # wgs vcf
parser.add_option("--anno_in", dest="anno_in")    # variant annotation

parser.add_option("--cadd_in", dest="cadd_in")    # CADD.tsv
parser.add_option("--cadd_indel_in", dest="cadd_indel_in")
parser.add_option("--dann_in", dest="dann_in")          # DANN core
parser.add_option("--chromHMM_in", dest="chromHMM_in")  # chromHMM
parser.add_option("--phylop_in", dest="phylop_in")      # phylop

# output
parser.add_option("--score_out", dest="score_out")  # Deleteriousness

(options, args) = parser.parse_args()

if os.path.exists(options.wgs_in) and os.path.exists(options.wgs_in+".tbi"):
    wgsTabix = pysam.Tabixfile(options.wgs_in,'r')                  # GTEx VCF file
if os.path.exists(options.af_in) and os.path.exists(options.af_in+".tbi"):
    afTabix = pysam.Tabixfile(options.af_in,'r')                    # GTEx AF file
if os.path.exists(options.anno_in) and os.path.exists(options.anno_in+".tbi"):
    annoTabix = pysam.Tabixfile(options.anno_in,'r')                    # variant annotation
if os.path.exists(options.cadd_in) and os.path.exists(options.cadd_in+".tbi"):
    cadd_snp_Tabix = pysam.Tabixfile(options.cadd_in,'r')                    # CADD.tsv
if os.path.exists(options.cadd_indel_in) and os.path.exists(options.cadd_indel_in+".tbi"):
    cadd_indel_Tabix = pysam.Tabixfile(options.cadd_indel_in,'r')
if os.path.exists(options.dann_in) and os.path.exists(options.dann_in+".tbi"):
    dannTabix = pysam.Tabixfile(options.dann_in,'r')                    # DANN score file
if os.path.exists(options.chromHMM_in) and os.path.exists(options.chromHMM_in+".tbi"):
    chromHMMTabix = pysam.Tabixfile(options.chromHMM_in,'r')            # chromHMM
if os.path.exists(options.phylop_in) and os.path.exists(options.phylop_in+".tbi"):
    phylopTabix = pysam.Tabixfile(options.phylop_in,'r')                # phylop



WGS_header = []
with gzip.open(options.wgs_in,"r") as f:
    for line in f:
        eachcol = line.rstrip().split('\t')
        if eachcol[0].startswith("#CHROM") == 1:
            for colnames in eachcol:
                WGS_header.append(colnames)
            break


dic_wgs = {}
for col_vcfs in WGS_header:
    dic_wgs[col_vcfs] = WGS_header.index(col_vcfs)



dic_chromHMM = {"1_Active_Promoter": 1, "2_Weak_Promoter": 2, "3_Poised_Promoter": 3, "4_Strong_Enhancer": 4, \
                "5_Strong_Enhancer": 4, "6_Weak_Enhancer": 5, "7_Weak_Enhancer": 5, "8_Insulator": 6, \
                "9_Txn_Transition": 7, "10_Txn_Elongation": 7, "11_Weak_Txn": 8, "12_Repressed": 9, \
                "13_Heterochrom/lo": 10, "14_Repetitive/CNV": 10, "15_Repetitive/CNV": 10}

dic_segway = {"C0": 1, "C1": 1, "D": 2, "E/GM": 3, "F0": 4, "F1": 4, "GE0": 5, "GE1": 5, "GE2": 5, "GM0": 6, "GM1": 6, \
              "GS": 7, "H3K9me1": 8, "L0": 9, "L1": 9, "R0": 10, "R1": 10, "R2": 10, "R3": 10, "R4": 10, "R5": 10, \
              "TF0": 11, "TF1": 11, "TF2": 11, "TSS": 12}

# Complete header from CADD data
cadd_header = []
with gzip.open(options.cadd_in,"r") as f:
    for line in f:
        eachcol = line.rstrip().split('\t')
        if eachcol[0].startswith("#Chrom") == 1:
            for colnames in eachcol:
                cadd_header.append(colnames)
            break

# from GC to Segway is from CADD file
score_header = ["Ensembl_id","Chrom","Pos","nvar","ref","alt","rv_allele","ind_ref_major","DistTSS","AF","AF_1KGeur","vcfGT","indID", \
		'anno', 'variant_class', \
		'Type', \
		'GC','CpG','priPhCons','mamPhCons','verPhCons','priPhyloP', \
                'mamPhyloP','verPhyloP','GerpN','GerpS','dnaHelT','dnaMGW','dnaProT','dnaRoll','fitCons', \
                'cHmmTssA','cHmmTssAFlnk','cHmmTxFlnk','cHmmTx','cHmmTxWk','cHmmEnhG','cHmmEnh','cHmmZnfRpts','cHmmHet', \
                'cHmmTssBiv','cHmmBivFlnk','cHmmEnhBiv','cHmmReprPC','cHmmReprPCWk','cHmmQuies','EncH3K27Ac','EncH3K4Me1', \
                'EncH3K4Me3','EncNucleo','EncOCCombPVal','EncOCDNasePVal','EncOCFairePVal','EncOCpolIIPVal','EncOCctcfPVal', \
                'EncOCmycPVal','EncOCDNaseSig','EncOCFaireSig','EncOCpolIISig','EncOCctcfSig','EncOCmycSig','TFBS','TFBSPeaks', \
                'TFBSPeaksMax','PHRED','Segway', \
		 "chromHMM",   "phylop",   "DANN"]
cadd_extract = ['Type', 'GC','CpG','priPhCons','mamPhCons','verPhCons','priPhyloP', \
                'mamPhyloP','verPhyloP','GerpN','GerpS','dnaHelT','dnaMGW','dnaProT','dnaRoll','fitCons', \
                'cHmmTssA','cHmmTssAFlnk','cHmmTxFlnk','cHmmTx','cHmmTxWk','cHmmEnhG','cHmmEnh','cHmmZnfRpts','cHmmHet', \
                'cHmmTssBiv','cHmmBivFlnk','cHmmEnhBiv','cHmmReprPC','cHmmReprPCWk','cHmmQuies','EncH3K27Ac','EncH3K4Me1', \
                'EncH3K4Me3','EncNucleo','EncOCCombPVal','EncOCDNasePVal','EncOCFairePVal','EncOCpolIIPVal','EncOCctcfPVal', \
                'EncOCmycPVal','EncOCDNaseSig','EncOCFaireSig','EncOCpolIISig','EncOCctcfSig','EncOCmycSig','TFBS','TFBSPeaks', \
                'TFBSPeaksMax','PHRED','Segway']



# gtex VEP format

##VEP=v83 cache=/humgen/atgu1/fs03/konradk/vep//homo_sapiens/83_GRCh37 db=. polyphen=2.2.2 sift=sift5.2.2 COSMIC=71 ESP=20141103 gencode=GENCODE 19 HGMD-PUBLIC=20152 genebuild=2011-04 regbuild=13 assembly=GRCh37.p13 dbSNP=144 ClinVar=201507
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info">

# udn VEP format

##VEP="v91" time="2018-03-19 10:24:27" cache="/users/xli6/data/xin/tools/ensembl-vep/vep_cache/homo_sapiens/91_GRCh37" ensembl=91.18ee742 ensembl-variation=91.c78d8b4 ensembl-funcgen=91.4681d69 ensembl-io=91.923d668 1000genomes="phase3" COSMIC="81" ClinVar="201706" ESP="20141103" HGMD-PUBLIC="20164" assembly="GRCh37.p13" dbSNP="150" gencode="GENCODE 19" genebuild="2011-04" gnomAD="170228" polyphen="2.2.2" regbuild="1.0" sift="sift5.2.2"
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE">


VEP_header = []
with gzip.open(options.anno_in,"r") as f:
    for line in f:
	match = re.search("^##INFO=<ID=CSQ,Number=\.,Type=String,Description=\"Consequence annotations from Ensembl VEP\. Format: ([\w_\|]*)\">", line)
	if match:
	    print match.group(1)
	    break

if not match:
    print "cannot find VEP format line"
    exit()

eachcol = match.group(1).rstrip().split('|')
for colnames in eachcol:
    VEP_header.append(colnames)

dic_vep = {}
for col_veps in VEP_header:
    dic_vep[col_veps] = VEP_header.index(col_veps)



# idx_labels = []
# for features in score_header[5:-3]:
#     if features == "DistTSS": continue # Distance to TSS is going to be added at the final step
#     else: idx_labels.extend([CADDtsv_header.index(features)])

score_out = open(options.score_out,'w');  score_out.write("\t".join(score_header)+"\n")

# read lines from text files for reading targeted regions
# chr | pos | ref | rv_allele | gencode_idx | ensembl_id | gene name | ind_major_ref | count_rv | distTSS
# ensembl_id | chr | pos
for target_region in sys.stdin:
    fields_site = target_region.rstrip().split('\t')

    dic_header = {}

    # (chrom, start, ref, alt, r_allele, gene_name, major==ref?, # var, DistTSS, AF, AF_1KG)
    chrom = int(fields_site[0])
    pos = int(fields_site[1])
    dic_header["Chrom"] = chrom
    dic_header["Pos"] = pos
    ref = fields_site[2]
    alt = fields_site[3]
    rv_allele = fields_site[4]
    dic_header["ref"] = ref
    dic_header["alt"] = alt
    dic_header["rv_allele"] = rv_allele
    dic_header["Ensembl_id"] = fields_site[5]
    ind_ref_major = int(fields_site[6]) # indicator of reference allele == major
    dic_header["ind_ref_major"] = ind_ref_major
    nvar = int(fields_site[7])
    dic_header["nvar"] = nvar
    DistTSS = int(fields_site[8])
    dic_header["DistTSS"] = DistTSS
    dic_header["AF"] = float(fields_site[9])
    dic_header["AF_1KGeur"] = float(fields_site[10])
    dic_header["vcfGT"] = fields_site[11]
    dic_header["indID"] = fields_site[12]

    count_feature = 13

    # background should be also n_var*background
    # out_score = [str(gene_name), str(idx_gencode), str(chrom), str(pos), str(n_var)]
    # print dic_header

    if (len(ref) == 1) and (len(rv_allele) == 1):
	caddTabix = cadd_snp_Tabix
    else:
	caddTabix = cadd_indel_Tabix

    # remove this debug portion to speed up    
    for poss_var in wgsTabix.fetch(chrom,pos-1,pos):
        fields_WGSvcf = poss_var.rstrip().split('\t')

	if pos != int(fields_WGSvcf[1]):
	    continue
	
        _ref = fields_WGSvcf[3]
	if ref != _ref :
	    continue

        set_allele = [fields_WGSvcf[3]]   # ref
        set_allele.extend(re.split(',',fields_WGSvcf[4])) # alt

	if (ref not in set_allele) or (rv_allele not in set_allele):
	    continue

        for target_site in afTabix.fetch(chrom,pos-1,pos):
            fields_af = target_site.rstrip().split('\t')
	    if int(fields_af[1]) != pos :
		continue
	    allele_af = fields_af[4].split(':')
	    if allele_af[0] != ref :
		continue
            minor_allele = []
	    all_allele = {}
	    ref_af = 0
	    max_af = 0
            for gtex_allele in fields_af[4:]:
                allele_af = gtex_allele.split(':')
                if  len(allele_af[0]) > 0:
		    if allele_af[0] == ref:
		    	ref_af = float(allele_af[1])
		    if float(allele_af[1]) > max_af:
			max_af = float(allele_af[1])
		    if allele_af[1] > 0:
			all_allele[allele_af[0]] = float(allele_af[1])
	    if ref_af == max_af:
		_ind_ref_major = 1
	    else:
		_ind_ref_major = 0
	    for one_allele in all_allele.keys():
		if all_allele[one_allele] < max_af:
		    minor_allele.extend([one_allele])
	

	if rv_allele not in minor_allele:
	    print "error extracting rare variants"
	    print dic_header
	    print target_site

	if _ind_ref_major != ind_ref_major:
	    print "discrepancy in indicating major allele, possible tie"
	    print dic_header
	    print target_site

	tempFORMAT = re.split(':', fields_WGSvcf[int(dic_wgs[options.ind_id])])
        gt_ind = tempFORMAT[0]
        # individual-specific allele
        gt_allele = re.split('[|/]',gt_ind)
        count_var = 0
	for gt in gt_allele[0:2]:
	    if gt == '.': continue
            elif int(gt) == set_allele.index(rv_allele):
		count_var += 1
	if count_var == 0:
	    continue #ambiguous genotype rows, two rows with same chrom, pos, {rv_allele, ...}, but different gentoypes for this individual
	if nvar != count_var:
	    print "error specifying allele counts"
	    print gt_allele, rv_allele, ind_ref_major, set_allele.index(rv_allele), nvar, count_var
	    print chrom, pos, fields_af
    """
    # remove this debug portion, to speed up
    list_tss = open('/users/xli6/data/xin/gtex/RIVER_input/reference/gencode.v19.genes.v7.patched_contigs.coding_lincRNA_autosome.bed','r')
    for tss_pos in list_tss.readlines():
        tss_info = tss_pos.rstrip().split('\t')
        if tss_info[0] == dic_header["Ensembl_id"]:
            tss1 = tss_info[2]
	    tss2 = tss_info[3]
	    strand = tss_info[4]
	
    if strand == '+':
    	_DistTSS = abs(int(tss1)-pos)
    else:
	_DistTSS = abs(int(tss2)-pos)
    if _DistTSS != DistTSS:
	print "error in calculating the TSS"
    """
        

    dic_header["anno"] = "undefined"
    dic_header["variant_class"] = "undefined"
    for feature in cadd_extract:
	dic_header[feature] = "NaN"
    dic_header["chromHMM"] = "NaN"
    dic_header["phylop"] = "NaN"
    dic_header["DANN"] = "NaN"



    if rv_allele != ref:
	list_anno = [];
        for anno_site in annoTabix.fetch(chrom,pos-1,pos):
            fields_anno = anno_site.rstrip().split('\t')

	    if pos != int(fields_anno[1]):
		continue
	    if ref != fields_anno[3]:
		continue
	    set_allele = [fields_anno[3]]   # ref
	    set_allele.extend(re.split(',',fields_anno[4])) # alt
	    if rv_allele not in set_allele:
		continue


	    # INFO fields are divided by ;
            # within CSQ each gene/allele is separated by ,
	    # then by | and by &
	    temp_field0 = re.split(';', fields_anno[7])
	    match = False
	    for temp_info in temp_field0:
		match = re.search("^CSQ=(.*)", temp_info)
		if match:
		    temp_field1 = match.group(1)
		    break
	    if not match:
		print "error in VEP annotation file, no CSQ information"
            temp_field2 = re.split('[,]',temp_field1) # per gene
	    # print temp_field2, rv_allele	    
            for fields_bygene in temp_field2[0:]:
                temp_field3 = re.split('[|]',fields_bygene)
		_allele = temp_field3[int(dic_vep['Allele'])]
		_gene = temp_field3[int(dic_vep['Gene'])]
		_allele_num = int(temp_field3[int(dic_vep['ALLELE_NUM'])])
		_variant_class = temp_field3[int(dic_vep['VARIANT_CLASS'])]
		_consequence = temp_field3[int(dic_vep['Consequence'])]
		vep_allele = _allele
		if rv_allele == set_allele[_allele_num]:
		    dic_header['variant_class'] = _variant_class
		    if _gene == re.split('[.]',dic_header["Ensembl_id"])[0]: # same gene
			list_anno.extend([_consequence])
		else:
		    continue
                    if vep_allele == "-":
                        if ref.startswith(rv_allele):
                            if temp_field3[1] not in list_anno:
                                list_anno.extend([_consequence])
                    if (vep_allele == rv_allele) or ((ref + vep_allele) == rv_allele):
                        if temp_field3[1] not in list_anno:
                            list_anno.extend([consequence])
		    
	
	count_feature += 1
	if len(list_anno) > 0:
	    dic_header["anno"] = "&".join(list_anno);
	else:
	    dic_header["anno"] = "undefined"

        for poss_site in caddTabix.fetch(chrom,pos-1,pos):
            #if poss_var.startswith('#'): continue
            fields_cadd = poss_site.rstrip().split('\t')
	    if pos != int(fields_cadd[1]):
		continue

            if (ref == fields_cadd[cadd_header.index("Ref")]) and (rv_allele == fields_cadd[cadd_header.index("Alt")]):   # match ref and alt allele
                for feature in cadd_extract:
                    if feature == "DistTSS": continue
                    elif feature == "PHRED":
                        dic_header[feature] = float(fields_cadd[cadd_header.index(feature)]); count_feature += 1;
                    elif fields_cadd[cadd_header.index(feature)] == "NA":
                        dic_header[feature] = "NaN"
                    elif feature == "Segway":
                        # dic_header[feature] = dic_segway[str(fields_cadd[cadd_header.index(feature)])]
			dic_header[feature] = str(fields_cadd[cadd_header.index(feature)]); 
			count_feature += 1
		    elif feature == "Type":
			dic_header[feature] = str(fields_cadd[cadd_header.index(feature)]); count_feature += 1
                    else:
                        dic_header[feature] = float(fields_cadd[cadd_header.index(feature)]); count_feature += 1
                break
            else: continue
        # DANN
        dic_header["DANN"] = "NaN"
        for dann_site in dannTabix.fetch(chrom,pos-1,pos):
            fields_dann = dann_site.rstrip().split('\t')
	    if pos != int(fields_dann[1]):
		continue
            if (fields_dann[2] == ref) and (fields_dann[3] == rv_allele):
                dic_header["DANN"] = float(fields_dann[4]); count_feature += 1

    elif rv_allele == ref: # ref is the rare allele
        dic_header["anno"] = "undefined"
	# in those cases rv_allele == ref, the effect is not defined for most annotation resources


        for poss_site in caddTabix.fetch(chrom,pos-1,pos):
            #if poss_var.startswith('#'): continue
            fields_cadd = poss_site.rstrip().split('\t')
	    if pos != int(fields_cadd[1]):
		continue
            if ref == fields_cadd[cadd_header.index("Ref")]:
                for feature in cadd_extract:
                    if feature == "DistTSS": continue
                    elif feature == "PHRED":
                        dic_header[feature] = "NaN"
                    elif fields_cadd[cadd_header.index(feature)] == "NA":
                        dic_header[feature] = "NaN"
                    elif feature == "Segway":
                        # dic_header[feature] = dic_segway[str(fields_cadd[cadd_header.index(feature)])]
			dic_header[feature] = str(fields_cadd[cadd_header.index(feature)])
			count_feature += 1
		    elif feature == "Type":
			dic_header[feature] = str(fields_cadd[cadd_header.index(feature)]); count_feature += 1
                    else:
                        dic_header[feature] = float(fields_cadd[cadd_header.index(feature)]); count_feature += 1
                break
        # DANN
        dic_header["DANN"] = "NaN"

    # chromHMM
    dic_header["chromHMM"] = "NaN"
    for poss_site in chromHMMTabix.fetch('chr'+str(chrom),pos-1,pos):
        fields_chromHMM = poss_site.rstrip().split('\t')
        # dic_header["chromHMM"] = dic_chromHMM[str(fields_chromHMM[3])]
	dic_header["chromHMM"] = str(fields_chromHMM[3])
	count_feature += 1

    # phylop
    dic_header["phylop"] = "NaN"
    for poss_site in phylopTabix.fetch('chr'+str(chrom),pos-1,pos):
        fields_phylop = poss_site.rstrip().split('\t')
        dic_header["phylop"] = float(fields_phylop[3]); count_feature += 1

    if count_feature > 1:
        out_score = []
        for feature in score_header:
            out_score.extend([str(dic_header[feature])])
        score_out.write("\t".join(out_score)+"\n")
	# print out_score



#print("--- %s hours ---" % ((time.time() - start_time)/3600))

#f = open('${RAREVARDIR}/RIVER/data/score/done.combined.txt', 'a+')
#f.write("%s\n" % (str(options.num)))
#f.close()







