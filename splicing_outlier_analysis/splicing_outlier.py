
import numpy as np
import sys
import glob, os
from scipy.stats import chisqprob, chisquare, chi2_contingency
from scipy import stats
from scipy.stats.distributions import chi2
import time
import argparse
from collections import defaultdict
from operator import truediv

junc_path='/srv/scratch/restricted/rare_diseases/data/splicing/juncfiles/all_filteredjunc/'

##-------- FUNCTIONS
## Get the junction file list
def make_juncfile_lists(input_text_file, rundir):   
	sys.stderr.write("\nGet list of files\n")
	list_of_files = []
	path_to_input_text_file="%s/%s"%(rundir,input_text_file)
	for junc in open(path_to_input_text_file):
		junc = junc.strip()
		if os.path.exists(rundir):
			list_of_files.append(junc)
		else: 
			sys.stderr.write("%s does not exist... check your junction files.\n"%junc)
			exit(0)
	return list_of_files



## This function aims at listing and grouping all junctions per chromosome from a list of files together and get the read coverage for each sample
def make_dictionary_junc(list_files):
	junctions= defaultdict(lambda: defaultdict(dict))
	for file_number, path_to_jxn_file in enumerate(list_files):
		sys.stderr.write("Processing junction file {0} \n".format(file_number))
		with open(path_to_jxn_file) as jxn_file:
			for junc_number, line in enumerate(jxn_file):
				#sys.stderr.write("Processing junction  {0} \n".format(junc_number))
				chrom, start, end, strand, intron_motif, annot, uniq, multi, max_overhang=line.rstrip().split('\t')
				junc=(start,end)
				#if chrom not in junctions:
				#	junctions[chrom]=[]
				if junc not in junctions[chrom]:
					junctions[chrom][junc] = len(list_files) * [0]
				junctions[chrom][junc][file_number]=int(uniq)
	return junctions



#This function get junctions from annotation file
def make_refjunc_dict(annot):
	dict_ref_junctions=defaultdict(list)
	with open(annot, 'r') as ref_junc:
		for line in ref_junc:
			chrom, start, end, gene=line.rstrip().split('\t')
			junction=(start,end)
			dict_ref_junctions[chrom].append(junction)
	return dict_ref_junctions


#This function tests for annotation status in a list of junctions
def is_junc_annot(junctions_list,dict_ref_junctions, chrom):
	statuses=[]
	for index, junction in enumerate(junctions_list):
		if junction in dict_ref_junctions[chrom]:
			statuses.append('A')
		else:statuses.append('NA')
	return statuses



# calculate ratio
def get_junc_ratio(junc_cov_list,sum_juncs_cov_list):
	ratio=map(truediv, junc_cov_list, sum_juncs_cov_list)
	return ratio


# this fuction calculates percentiles of interest for jucntion ratios
def get_percentiles(junc_ratio_values,list_of_percentiles=[1,2,5,10,90,95,98,99]):
	percentiles=np.percentile(junc_ratio_values, q=list_of_percentiles)
	return percentiles

# this function returns outlier indexes for each junction ratio
def get_outlier_indexes(junc_ratio_values,junc_percentiles_values):
	outliers={'1':[],'2':[],'5':[],'10':[],'90':[],'95':[],'98':[], '99':[]}
	for indiv, ratio in enumerate(junc_ratio_values):
		if ratio<=junc_percentiles_values[0]:
			outliers['1'].append(indiv)
		if ratio<=junc_percentiles_values[1]:
			outliers['2'].append(indiv)
		if ratio<=junc_percentiles_values[2]:
			outliers['5'].append(indiv)
		if ratio<=junc_percentiles_values[3]:
			outliers['10'].append(indiv)
		if ratio>=junc_percentiles_values[4]:
			outliers['90'].append(indiv)
		if ratio>=junc_percentiles_values[5]:
			outliers['95'].append(indiv)
		if ratio>=junc_percentiles_values[6]:
			outliers['98'].append(indiv)
		if ratio>=junc_percentiles_values[7]:
			outliers['99'].append(indiv)
	return outliers

# this fuction calculates zscores for junction ratios
def get_zscores(junc_ratio_values):
	zscores=stats.zscore(np.array(junc_ratio_values))
	return zscores


# This function returns index of samples that have |zscore| >=2 
def get_outlier_indexes_zscore(junc_zscore_values,threshold=2):
	outliers_zscore=[]
	for indiv, ratio in enumerate(junc_zscore_values):
		if abs(ratio)>=threshold:
			outliers_zscore.append(indiv)
	return outliers_zscore


#This function makes a dictionnary with samples as values and affected status as key
def make_affectstatus_dict(meta_file):
	affectstatus_dict={'Case':[],'Control':[]}
	with open(meta_file, 'r') as meta:
		#next(meta)
		for line in meta:
			sample_id,affected_status=line.rstrip().split('\t')
			if affected_status == "Case":
				affectstatus_dict['Case'].append(sample_id)
			if affected_status == "Control":
				affectstatus_dict['Control'].append(sample_id)
	return affectstatus_dict

# This function makes a list of samples ordered as in the list_of_files processed
def make_sample_list(list_of_files):
	file_temp=[i.split(junc_path, 1)[1] for i in list_of_files]
	samples=[i.split('.', 1)[0] for i in file_temp]
	return samples

# this function counts the number of outlier in cases and controls
def outlier_status(samples_list,affected_status_dict):
	status_dic={'Case':0,'Control':0}
	for sample in samples_list:
		if sample in affected_status_dict['Case']:
			status_dic['Case']+=1
		elif sample in affected_status_dict['Control']:
			status_dic['Control']+=1
	return status_dic.values()


# get sample id from indexes
# get affected status from sample ID
# count number of cases and controls
# return affected status of outliers for each percentile

#for each key of sample_ids_dic , get affected status for each element of the value
def get_outlier_affected_status(junc_outlier_values_dict,samples,affected_status_dict):
	sample_ids_dic={k: list(np.array(samples)[v]) for k, v in  junc_outlier_values_dict.items()}
	affected_status_dic={k:outlier_status(v,affected_status_dict) for k,v in sample_ids_dic.items()}
	return  affected_status_dic

# this function returns affected status
def get_outlier_zscore_affected_status(outlier_index_list,samples,affected_status_dict):
# for each list get the sample name
	sample_ids=[np.array(samples)[index] for index in  outlier_index_list]
# get the affected status and report number of case vs controls
	affected_status=outlier_status(sample_ids,affected_status_dict)
	return  affected_status

# this function is calculating the coverage prop for each junction for a single donor or acceptor
def process_junctions(annot,junctions_dict,dict_ref_junctions,direction, rundir, outPrefix,list_of_files, samples,affected_status_dict):
	outFile = "%s/%s_junc_outliers.txt"%(rundir,outPrefix)
	outratio="%s/%s_junc_ratios.txt"%(rundir,outPrefix)
	with open(annot, 'r') as ref_junc, open(outFile, 'w') as OUT, open(outratio, 'w') as OUTratio:
		#header= ['gene', 'chr', 'junction', 'annotation_status','Cluster', 'OutlierCase', 'OutlierControl']+samples#+[(str(e)) for e in range(len(list_of_files))]
		header= [ 'chr', 'junction_start', 'junction_end','gene', 'annotation_status','Cluster', 'OutlierCase', 'OutlierControl']+samples
		OUT.write('\t'.join(header) + '\n')
		header_ratio=['chr', 'junction_start', 'junction_end','gene', 'annotation_status','Cluster']+samples
		OUTratio.write('\t'.join(header_ratio) + '\n')
		for ref_junc_num, line in enumerate(ref_junc):
			chrom, start, end, gene=line.rstrip().split('\t')
			if direction == 'donor':
				valid_jxns= filter(lambda jxn: jxn[0] == start, junctions_dict[chrom].keys())
			if direction == 'acceptor':
				valid_jxns= filter(lambda jxn: jxn[1] == end, junctions_dict[chrom].keys())
			if len(valid_jxns)==1:
				continue
			statuses=is_junc_annot(valid_jxns,dict_ref_junctions, chrom)
			test_dic={}
			for jxn in valid_jxns:
				test_dic[jxn] = junctions_dict[chrom].get(jxn, len(list_of_files)*[0])
				#test_dic_ctrl[jxn] = junc_control[chrom].get(jxn, len(list_of_controls_files)*[0])
			# add one to every count before computing ratio
			#test_dic={k: np.array(v)+1 for  k, v in test_dic.items()}
			# get sum of counts over junctions
			junc_sum=[sum(x)for x in zip(*test_dic.values())]
			# filter samples always at 0
			junc_sum_filt= [x if x != 0 else np.nan for x in junc_sum]
			# calculate the ratio of counts for each junction and each individual
			junc_ratio_dic = {k: get_junc_ratio(v,junc_sum_filt) for k, v in test_dic.items()}
			# get percentiles for each junction of teh dictionary
			#junc_percentiles_dic={k: get_percentiles(v) for k, v in junc_ratio_dic.items()}
			# get outlier indexes
			#junc_outlier_dic={k: get_outlier_indexes(v, junc_percentiles_dic[k]) for k, v in junc_ratio_dic.items()}
			# get number of case and controls for outliers
			#junc_affected_status_dic={k: get_outlier_affected_status(v, samples,affected_status_dict) for k, v in junc_outlier_dic.items()}
			# get zcores for ratios
			junc_zscore_dic={k: get_zscores(v) for k, v in junc_ratio_dic.items()}
			# get zscore outlier
			junc_outlier_zscore_dic={k: get_outlier_indexes_zscore(v) for k, v in junc_zscore_dic.items()}
			# count affected status for zscore outliers
			junc_affected_status_zscore={k: get_outlier_zscore_affected_status(v,samples,affected_status_dict) for k, v in junc_outlier_zscore_dic.items()}
			# report results for junction set
			for key_number, key in enumerate(junc_ratio_dic):
				values=[chrom]+[i for i in key]+[gene, statuses[key_number], ref_junc_num]+[outlier_number for outlier_number in junc_affected_status_zscore[key]] +junc_zscore_dic[key].tolist()
				values_ratio=[chrom]+[i for i in key]+[gene, statuses[key_number], ref_junc_num]+junc_ratio_dic[key]
				#values=[gene, chrom,key, statuses[key_number], ref_junc_num]+[outlier_number for outlier_number in junc_affected_status_zscore[key]] +junc_zscore_dic[key].tolist()# [value for value in junc_ratio_dic[key]]
				OUT.write('\t'.join([str(s) for s in values])+ '\n')
				OUTratio.write('\t'.join([str(s) for s in values_ratio])+ '\n')



## This function is filtering output for unique lines after both donors and acceptors were tested
def filter_outfile(outPrefix, rundir):
	inFile = "%s/%s_junc_outliers.txt"%(rundir,outPrefix)
	outFile="%s/%s_junc_outliers_filtered.txt"%(rundir,outPrefix)
	inFile_ratio="%s/%s_junc_ratios.txt"%(rundir,outPrefix)
	outFile_ratio="%s/%s_junc_ratios_filtered.txt"%(rundir,outPrefix)
	lines_seen = set()
	lines_seen_ratio=set() # holds lines already seen
	outfile = open(outFile, "w")
	outfile_ratio=open(outFile_ratio, 'w')
	for line in open(inFile, "r"):
		if line not in lines_seen: # not a duplicate
			outfile.write(line)
			lines_seen.add(line)
	outfile.close()
	for line in open(inFile_ratio, "r"):
		if line not in lines_seen_ratio: # not a duplicate
			outfile.write(line)
			lines_seen_ratio.add(line)
	outfile_ratio.close()





def main():
	USAGE = """
		Takes file with list of control junction files, and a file with list of case junction files. 
		Test junctions from annotation and look for different slice events between controls and case.
    	"""

	parser= argparse.ArgumentParser(description = USAGE)
	parser.add_argument("--juncfiles", dest="juncfiles",required=True,help="text file with all junction files to be processed")
	parser.add_argument("--outprefix", dest="outprefix", default = 'RD', help="output prefix (default RD)")
	parser.add_argument("--rundir", dest="rundir", default='./',help="write to directory (default ./)")
	parser.add_argument("--annot", dest="annot_junc", default='/srv/scratch/restricted/rare_diseases/data/splicing/juncfiles/gencodev19_intronsstartplus1_proteincoding.genenames_uniqjunc.tsv', help="file containing known annotated junctions (default /srv/scratch/restricted/rare_diseases/data/splicing/juncfiles/gencodev19_intronsstartplus1_proteincoding.genenames_uniqjunc.tsv)")
	parser.add_argument("--meta_file", dest="meta_file", default='/srv/scratch/restricted/rare_diseases/data/metadata/affected_status_freeze.txt',help="file containing sample_id and affected status")
	
	args= parser.parse_args()


	list_junctions_file=args.juncfiles
	outPrefix = args.outprefix
	rundir = args.rundir
	annot=args.annot_junc
	meta_file=args.meta_file

	sys.stderr.write(time.ctime())

	print 'Create list of files'
	list_of_files=make_juncfile_lists(list_junctions_file, rundir)

	print 'Create dictionaries'
	junctions_dict=make_dictionary_junc(list_of_files)

	print 'Make annotated junctions dictionary'
	dict_ref_junctions=make_refjunc_dict(annot)

	print 'Make dict of case and control samples'
	affected_status_dict=make_affectstatus_dict(meta_file)
	
	print 'Get sample names'
	samples=make_sample_list(list_of_files)

	print 'Get ratios for donors'
	process_junctions(annot,junctions_dict,dict_ref_junctions, 'donor', rundir, outPrefix,list_of_files,samples,affected_status_dict)
	


	print 'Get ratios for acceptors'
	process_junctions(annot,junctions_dict,dict_ref_junctions, 'acceptor', rundir, outPrefix,list_of_files,samples,affected_status_dict)

	print 'Filter results for unique lines'
	filter_outfile(outPrefix, rundir)
	
	sys.stderr.write(time.ctime())




if __name__ == "__main__":
	main()


