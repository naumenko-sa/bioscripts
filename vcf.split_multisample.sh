#!/bin/bash

# splits multisample vcfs in the current directory
# https://www.biostars.org/p/138694/

for file in *.vcf*; 
do
    for sample in `bcftools query -l $file`; 
    do
	bcftools view -c1 -Oz -s $sample -o ${file/.vcf*/.$sample.vcf.gz} $file
    done
done
