#!/bin/env python

# $1 = file.bam.coverage = product of ~/bioscripts/bam.coverage.sh with -d on

import sys

cov={}

levels = [0,1,10,20,30,50,100,200]
for i in levels:
    cov[i]=0

total_bases=0

with open(sys.argv[1],'rb') as f_coverage:
    for line in f_coverage:
	total_bases += 1
	buf = line.strip()
	fields = buf.split('\t')
	coverage = int(fields[5])

	for j in levels:
	    if (coverage >= j):
		cov[j] += 1

f_coverage.close()

for i in levels:
    print(str(i)+'\t'+str(cov[i])+'\t'+str(1.0*cov[i]/total_bases))
