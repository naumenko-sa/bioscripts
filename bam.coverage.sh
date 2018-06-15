#!/bin/bash
# $bam - input.bam
# $bed - filter.bed
# output - input.bam.coverage
# http://bedtools.readthedocs.io/en/latest/content/tools/coverage.html

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

#-d gives coverage of every nucleotide
#without sorted it gives alloc error
#for sorted the order of chromosomes in a bed and in a bam file should be the same
#https://github.com/arq5x/bedtools/issues/109
#bedtools index just chr name and length
#don't forget to sort bed : bedtools sort -faidx names.txt
bedtools coverage -d -sorted -a $bed -b $bam -g /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa.bedtoolsindex > $bam.dcoverage
