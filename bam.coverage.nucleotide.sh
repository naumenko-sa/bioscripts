#!/bin/bash
# calculates coverage of nucleotides in a bed file and a bam file
# to work with the particular gene, subset a bam file first:
# samtools view -b DMD.bam X > DMD.X.bam

# $bam - input.bam
# $bed - filter.bed
# output - input.bam.coverage

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

#uses bedtools from bcbio

bedtools coverage -a $bed -b $bam -d -sorted -g /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa > $bam.coverage
