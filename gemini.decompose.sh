#!/bin/bash

# looks like bcbio uses https://github.com/vcflib/vcflib#vcfallelicprimitives for multiallic SNPs decomposition, and finally gemini has problems with some fields
# VT should be better: https://github.com/atks/vt as gemini site suggests: 
# https://gemini.readthedocs.io/en/latest/content/preprocessing.html

# $1 = test.vcf

bname=`echo $1 | sed s/.vcf//`
cat $1 | sed s/"##FORMAT=<ID=AD,Number=."/"##FORMAT=<ID=AD,Number=R"/ > $bname.modified.vcf
vt decompose -s $bname.modified.vcf | vt normalize -r ~/work/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa - | vt uniq - > $bname.decomposed.vcf

