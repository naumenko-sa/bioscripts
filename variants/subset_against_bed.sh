#!/bin/bash

# $1 - vcf.gz
# $2 - bed
# #3 - suffix for the output vcf

bname=`basename $1 .vcf.gz`

if [ $# -lt 3 ];then
    suffix="in_panel"
else
    suffix=$3
fi

# extract header
gunzip -c $1 | grep "^#" > $bname.$suffix.vcf
bedtools intersect -a $1 -b $2 >> $bname.$suffix.vcf
bgzip $bname.$suffix.vcf
tabix $bname.$suffix.vcf.gz
