#!/bin/bash
# split multisample vcf
# https://www.biostars.org/p/138694/

for sample in `bcftools query -l $file`; 
do
    bcftools view -c1 -Oz -s $sample -o ${file/.vcf*/.$sample.vcf.gz} $file
done
