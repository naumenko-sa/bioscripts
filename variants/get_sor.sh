#!/bin/bash

# for a vcf sample, get SOR values for HET and HOM sites for SNVs
sample=`basename $1 .vcf.gz`

gunzip -c $1 | grep -v "^#" | grep PASS | grep "0/1"  | awk -F '\t' '{if(length($4) == length($5)) print $8}' | awk -F '=' '{print $NF}' > $sample.het.sor.tsv
gunzip -c $1 | grep -v "^#" | grep PASS | grep "1/1"  | awk -F '\t' '{if(length($4) == length($5)) print $8}' | awk -F '=' '{print $NF}' > $sample.hom.sor.tsv
