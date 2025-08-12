#!/bin/bash

# $1 - clinvar.PANELWES.missing.vcf.gz
# get MC stat
bname=`basename $1 .vcf.gz`
gunzip -c $1 | grep -v "^#" |  awk -F 'MC=' '{print $2}' | awk -F ';' '{print $1}' | sort  | uniq -c | awk '{print $2"\t"$1}' | sort -k2,2nr > $bname.MC_stat.tsv

