#!/bin/bash

# bcbio uses vt decompose but not normalize and uniq
# https://gemini.readthedocs.io/en/latest/content/preprocessing.html

# $1 = test.vcf.gz

#bgzip -d $1
bname=`basename $1 .vcf.gz`

#already fixed in the new version
#cat $bname.vcf | sed s/"##FORMAT=<ID=AD,Number=."/"##FORMAT=<ID=AD,Number=R"/ > $bname.modified.vcf

vt decompose -s $1 | vt normalize -r ~/work/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa -n - | vt uniq - -o $bname.decomposed.vcf.gz

tabix $bname.decomposed.vcf.gz
