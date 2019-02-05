def getTerm(stream):
	block = []
	for line in stream:
		if line.strip() == "[Term]" or line.strip() == "[Typedef]":
			break
		else:
			if line.strip() != "":
				block.append(line.strip())
	return block

def parseTagValue(term):
	data = {}
	for line in term:
		tag = line.split(': ',1)[0]
		value = line.split(': ',1)[1]
		if not data.has_key(tag):
			data[tag] = []
		data[tag].append(value)
	return data

oboFile = open('18_10_23_hp.obo','r')
#declare a blank dictionary
#keys are the HPO ids
terms = {}

#skip the file header lines
getTerm(oboFile)


#infinite loop to go through the obo file.
#Breaks when the term returned is empty, indicating end of file
while 1:
	#get the term using the two parsing functions
	term = parseTagValue(getTerm(oboFile))
	if len(term) != 0:
		termID = term['id'][0]
		#only add to the structure if the term has a is_a tag
		#the is_a value contain GOID and term definition
		#we only want the GOID
		if term.has_key('is_a'):
			termParents = [p.split()[0] for p in term['is_a']]
			if not terms.has_key(termID):
				#each goid will have two arrays of parents and children
				terms[termID] = {'p':[],'c':[], 'a':[]}
			#append parents of the current term
			terms[termID]['p'] = termParents
			#for every parent term, add this current term as children
			for termParent in termParents:
				if not terms.has_key(termParent):
					terms[termParent] = {'p':[],'c':[],'a':[]}
				terms[termParent]['c'].append(termID)
		if term.has_key('alt_id'):
			termAlt=[p.split()[0] for p in term['alt_id']]
			if not terms.has_key(termID):
				#each goid will have two arrays of parents and children
				terms[termID] = {'p':[],'c':[], 'a':[]}
			terms[termID]['a']=termAlt
	else:
		break

def getDescendents(goid):
	recursiveArray = [goid]
	if terms.has_key(goid):
		children = terms[goid]['c']
		if len(children) > 0:
			for child in children:
				recursiveArray.extend(getDescendents(child))

	return set(recursiveArray)

def getAncestors(goid):
	recursiveArray = [goid]
	if terms.has_key(goid):
		parents = terms[goid]['p']
		if len(parents) > 0:
			for parent in parents:
				recursiveArray.extend(getAncestors(parent))

	return set(recursiveArray)


with open('parent_terms_hpo_2018_10_23.txt', 'w') as f:
	for k, v in terms.iteritems():
		myline=[k]+ [','.join(str(p) for p in v['p'])]
		f.write('\t'.join([str(s) for s in myline])+ '\n')

with open('child_terms_hpo_2018_10_23.txt', 'w') as f:
	for k, v in terms.iteritems():
		myline=[k]+ [','.join(str(p) for p in v['c'])]
		f.write('\t'.join([str(s) for s in myline])+ '\n')


with open('altid_terms_hpo_2018_10_23.txt', 'w') as f:
	for k, v in terms.iteritems():
		myline=[k]+ [','.join(str(p) for p in v['a'])]
		f.write('\t'.join([str(s) for s in myline])+ '\n')

