#!/bin/bash

# $1 - vcf.gz
# $2 - bed

bname=`basename $1 .vcf.gz`

# extract header
gunzip -c $1 | grep "^#" > $bname.in_panel.vcf
bedtools intersect -a $1 -b $2 >> $bname.in_panel.vcf
bgzip $bname.in_panel.vcf
tabix $bname.in_panel.vcf.gz
