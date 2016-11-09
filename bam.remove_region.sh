#!/bin/bash

#removes reads associated with a specific region from a bam file

# $bam - input.bam
# $bed - filter.bed
# output - input.filtered.bam

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

module load bedtools

bedtools intersect -abam $bam -b $bed -v > `echo $bam | sed s/bam/filtered.bam/`
