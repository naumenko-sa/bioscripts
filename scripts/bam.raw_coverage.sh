#!/bin/bash
# $bam - input.bam
# $bed - filter.bed
# output - input.bam.coverage

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

if [ -z $bam ]
then
    bam=$1
fi

if [ -z $bed ]
then
    bed=$2
fi

bedtools multicov -bams $bam -bed $bed > $bam.raw_coverage
