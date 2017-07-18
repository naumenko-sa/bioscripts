#!/bin/bash
# get coverage for exons of a gene
# $sample - sample name, $sample.chrom.bam should be in the current directory
# $gene - gene, gene.bed should be in the current directory
# output - $sample.$gene.coverage
# first generate bam files for each chromosome:
# samtools view -bh DMD.bam X > DMD.X.bam
# samtools index DMD.X.bam
# otherwise (1chr for bed, 1chr for bam it is not working)
# -split for RNA-seq coverage

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

chrom=`head -n1 $bed | cut -f1`
sample=`echo $bam | awk -F "." '{print $1}'`

echo "Calculates coverage of" $gene.bed " in " $sample " for chr " $chr

echo -e "chrom\tstart\tend\texon\t$sample" > $bam.$bed.coverage
bedtools coverage -a $bed -b $bam -split -mean -sorted >> $bam.$bed.coverage
