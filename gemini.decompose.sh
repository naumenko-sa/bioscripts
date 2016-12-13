#!/bin/bash

# bcbio uses vt decompose but not normalize and uniq
# https://gemini.readthedocs.io/en/latest/content/preprocessing.html

# $1 = test.vcf.gz

bgzip -d $1
bname=`echo $1 | sed s/.vcf.gz//`
cat $bname.vcf | sed s/"##FORMAT=<ID=AD,Number=."/"##FORMAT=<ID=AD,Number=R"/ > $bname.modified.vcf
vt decompose -s $bname.modified.vcf | vt normalize -r ~/work/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa - | vt uniq - > $bname.decomposed.vcf
bgzip $bname.vcf
bgzip $bname.decomposed.vcf
tabix $bname.decomposed.vcf.gz
rm $bname.modified.vcf

