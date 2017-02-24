#!/bin/bash

#for CHEO project:
#create filtered vcf from vcf with all variants and the final report in csv format

#$1 - project name

# variants in project-ensemble.vcf.gz are not normalized (multiallelic snps are in the same record), 
# in final report they are
# use uniq to avoid problems

# no need to sort - it will change the order of SNPs
# normalized records should follow each other

cd $1
cat $1.csv  | sed 1d | awk -F ',' '{print $1}' | sed s/'"'//g | uniq | sed s/":"/"\t"/ | sed s/chr// > $1.index
bcftools view -R $1.index -O z -o $1.vcf.gz ${1}-ensemble.vcf.gz
rm $1.index
cd ..