#!/bin/bash
# $bam - input.bam
# $bed - panel.bed
# output - input.bam.coverage
# http://bedtools.readthedocs.io/en/latest/content/tools/coverage.html

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=2g,mem=2g

if [ -z $bam ]
then
    bam=$1
fi

if [ -z $bed ]
then
    bed=$2
fi

if [ -z $type ]
then
   type=$3
fi

#histogram
#bedtools coverage -hist -a $bam -b $bed > $bam.coverage

#-d gives coverage of every nucleotide
#without sorted it gives alloc error
#for sorted the order of chromosomes in a bed and in a bam file should be the same
#https://github.com/arq5x/bedtools/issues/109
#bedtools index just chr name and length
#don't forget to sort bed : bedtools sort -faidx names.txt

params=''

if [ "$type" == "rnaseq" ]
then
    params=' -split'
fi

echo "Start: " `date`

#bedtools sort -faidx /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa.bedtoolsindex -i $bed > $bed.faidxsorted.bed

bedtools coverage -d -sorted -a $bed -b $bam -g /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa.bedtoolsindex $params > $bam.dcoverage

bam.coverage.base_wise.stat.py $bam.dcoverage > $bam.coverage_stats

median_line=`cat $bam.dcoverage | wc -l`
median_line=$(($median_line/2))
cat $bam.dcoverage | awk '{print $6}' | sort -n | sed -n ${median_line}p > $bam.median

echo "END: " `date`
