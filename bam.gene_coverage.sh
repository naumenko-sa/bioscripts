#!/bin/bash
# get coverage for a exons in a gene
# $sample - sample name, $sample.chrom.bam should be in the current directory
# $gene - gene, gene.bed should be in the current directory
# output - $sample.$gene.coverage
# first generate bam files for each chromosome:
# samtools view -bh DMD.bam X > DMD.X.bam

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=15g,mem=15g

if [ -z $sample ]
then
    sample=$1
fi

if [ -z $gene ]
then
    gene=$2
fi

chrom=`cat $gene.bed | head -n1 | cut -f1`

echo $sample $gene

echo -e "chrom\tstart\tend\texon\t$sample" > $sample.$gene.coverage
bedtools coverage -a $gene.bed -b $sample.$chrom.bam -mean >> $sample.$gene.coverage
