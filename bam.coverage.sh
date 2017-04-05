#!/bin/bash
#calculates coverage of a bed file with mapped reads from bam file
#not very usable, too low level: outputs all reads
#use bamstats04 instead from jvarkit

# $bam - input.bam
# $bed - filter.bed
# output - input.bam.coverage

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

#module load bedtools
#uses bedtools from bcbio

if [ -z $bam ]
then
    bam=$1
fi

if [ -z $bed ]
then
    bed=$2
fi

bedtools coverage -hist -abam $bam -b $bed > $bam.coverage
