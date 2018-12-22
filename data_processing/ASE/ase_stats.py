import sys
import numpy as np
import pandas as pd
import scipy.stats
import glob

# usage: python ase_stats.py [infile] [outfile]
#modified by Nikki Teran from Nicole Ferraro Nov 30 2017
#files = glob.glob('../output/filtered_ase_files/*ase_out_filter.txt')

#From $RARE_DIS_DIR/ase/scripts/ase_model.py
### Read data
f = sys.argv[1]
out = sys.argv[2]

data = pd.read_csv(f, sep='\t')

### Empirically estimate prior beta distribution
alpha0, beta0, _, _ = scipy.stats.beta.fit(data['refRatio'].values, 1, 1, floc=0, fscale=1)

### Perform shrinkage
data['eb_ref_ratio'] = pd.Series((data['refCount'] + alpha0) / (data['rawDepth'] + alpha0 + beta0), index=data.index)
data['alpha1'] = pd.Series(data['refCount'] + alpha0, index=data.index)
data['beta1'] = pd.Series(data['rawDepth'] - data['refCount'] + beta0, index=data.index)

low, high = scipy.stats.beta.interval(.95, data['alpha1'], data['beta1'])

data['eb_low'] = pd.Series(low, index=data.index)
data['eb_high'] = pd.Series(high, index=data.index)

### Point estimate
pep1 = scipy.stats.beta.cdf(alpha0/(alpha0+beta0), data['alpha1'], data['beta1'])
pep2 = 1-scipy.stats.beta.cdf(alpha0/(alpha0+beta0), data['alpha1'], data['beta1'])
data['PEP'] = pd.Series(np.minimum(pep1, pep2), index=data.index)
sorted_keys = list(np.argsort(data['PEP'].values))

sorted_values = np.sort(data['PEP'])
qvalues = np.cumsum(sorted_values) / np.arange(1, len(sorted_values)+1)
values = np.ones(len(qvalues))

values[sorted_keys] = qvalues
data['qvalue'] = pd.Series(values, index=data.index)

data.to_csv(out, sep='\t', index=False)
