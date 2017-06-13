#!/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/anaconda/bin/python

import csv
import sys
from collections import defaultdict

fieldnames=['Position','Ref','Alt','Variation','Zygosity','Protein_change_ensembl','Gene','Conserved_in_29_mammals','Sift_score','Polyphen_score','Cadd_score']
#Frequency','Samples']

frequencies = defaultdict(list)
samples = defaultdict(list)
annotations = defaultdict(list)

with open(sys.argv[1],'r') as f_csv:
    reader = csv.DictReader(f_csv)
    for row in reader:
	superkey=row['Position']+'-'+row['Ref']+'-'+row['Alt']
	if superkey in frequencies:
	    frequencies[superkey] += 1
	    samples[superkey].append(row['Sample'])
	else:
	    frequencies[superkey] = 1
	    l = []
	    for key in fieldnames:
		l.append(row[key])
	    annotations[superkey] = ','.join(l)
	    ll = []
	    ll.append(row['Sample'])
	    samples[superkey] = ll

with open(sys.argv[2],'w') as f_out:
    f_out.write(','.join(fieldnames)+',Frequency,Samples')
    f_out.write('\n')
    for key in sorted(frequencies.keys()):
	f_out.write(annotations[key]+','+str(frequencies[key])+','+';'.join(samples[key]))
	f_out.write('\n')

f_out.close()
f_csv.close()
