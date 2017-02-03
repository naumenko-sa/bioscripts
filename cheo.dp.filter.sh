#!/bin/bash

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

#parameters:
#vcf
#DP
#filters out variants with DP < $DP
#uses gatk wrapper from bcbio

gatk -R /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa \
     -T VariantFiltration \
     -V $vcf \
     -o `echo $vcf | sed s/vcf/DPLT$DP.vcf/` \
     --filterExpression " DP < $DP" --filterName "DPLT$DP"