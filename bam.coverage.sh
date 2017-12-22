#!/bin/bash
# $bam - input.bam
# $bed - filter.bed
# output - input.bam.coverage

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=20g,mem=20g

if [ -z $bam ]
then
    bam=$1
fi

if [ -z $bed ]
then
    bed=$2
fi

#histogram
#bedtools coverage -hist -a $bam -b $bed > $bam.coverage

#coverage of every nucleotide
bedtools coverage -d -a $bed -b $bam > $bam.dcoverage
