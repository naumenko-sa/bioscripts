#!/bin/bash

# creates gemini database with vep refseq annotations, 
# dumps text file and variant_impacts for rare harmful variants

#PBS -l walltime=10:00:00,nodes=1:ppn=10
#PBS -joe .
#PBS -d .
#PBS -l vmem=50g,mem=50g


#vt rminfo ${family}-ensemble.vcf.gz -t CSQ -o ${family}.no_vep.vcf.gz

#gemini.decompose.sh ${family}.no_vep.vcf.gz

#gemini.vep.refseq.sh ${family}.no_vep.decomposed.vcf.gz

gemini.vep2gemini.sh ${family}.no_vep.decomposed.vepeffects_refseq.vcf.gz

gemini.gemini2txt.sh ${family}.no_vep.decomposed.vepeffects_refseq.db

gemini.variant_impacts.sh ${family}.no_vep.decomposed.vepeffects_refseq.db
