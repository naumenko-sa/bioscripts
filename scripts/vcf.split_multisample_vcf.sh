#!/bin/bash
# split multisample vcf
# https://www.biostars.org/p/138694/

bname=`basename $1 .vcf.gz`


for sample in `bcftools query -l $1`
do
    bcftools view -c1 -Oz -s $sample -o $bname.$sample.vcf.gz $1
done
