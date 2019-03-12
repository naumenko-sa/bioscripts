#!/bin/bash

# creates gemini database with vep annotations and all impacts
# it is for projects that were done with previous version - with one impact per gene
# dumps text file and variant_impacts for rare harmful variants

#PBS -l walltime=10:00:00,nodes=1:ppn=10
#PBS -joe .
#PBS -d .
#PBS -l vmem=50g,mem=50g

if [ -z $family ];
then
    family=$1
fi

vt rminfo ${family}-ensemble.vcf.gz -t CSQ -o ${family}.no_vep.vcf.gz

gemini.decompose.sh ${family}.no_vep.vcf.gz

gemini.vep.sh ${family}.no_vep.decomposed.vcf.gz

gemini.vep2gemini.sh ${family}.no_vep.decomposed.vepeffects.vcf.gz

gemini.gemini2txt.sh ${family}.no_vep.decomposed.vepeffects.db

gemini.variant_impacts.sh ${family}.no_vep.decomposed.vepeffects.db
