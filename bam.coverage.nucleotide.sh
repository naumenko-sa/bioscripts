#!/bin/bash
# calculates coverage of nucleotides in a bed file for all bam files in the current directory

# $bam - input.bam
# $bed - filter.bed
# output - filter.bed.coverage

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

#module load bedtools
#uses bedtools from bcbio

bedtools coverage -a $bed -b *.bam -d -sorted > $bed.coverage
