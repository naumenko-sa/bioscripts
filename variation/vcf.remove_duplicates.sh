#!/bin/bash

bname=`basename $1 .vcf.gz`

# bcftools removing alleles in muliallelic sites
#bcftools norm --rm-dup snps -Oz $1 > $bname.no_duplicates.vcf.gz

vt uniq $1 | bgzip > $bname.no_duplicates.vcf.gz
tabix $bname.no_duplicates.vcf.gz
tabix --csi $bname.no_duplicates.vcf.gz

