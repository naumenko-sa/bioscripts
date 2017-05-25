#!/bin/bash
# get coverage for a exons in a gene
# $bam - input.bam
# $bed - gene.bed
# output - input.bam.coverage
# first generate bam files for each chromosome:
# samtools view -bh DMD.bam X > DMD.X.bam

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

if [ -z $sample ]
then
    sample=$1
fi

if [ -z $bed ]
then
    bed=$2
fi

chrom=`cat $bed | head -n1 | cut -f1`

bedtools coverage -a $bed -b $sample.$chrom.bam -mean > $sample.$bed.coverage
