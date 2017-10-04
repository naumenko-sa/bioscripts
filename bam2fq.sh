#!/bin/bash

#PBS -l walltime=10:00:00,nodes=1:ppn=10
#PBS -joe .
#PBS -d .
#PBS -l vmem=20g

if [ -z $bam ]
then
    bam=$1
fi

if [ -z $sample ]
then
    sample=$2
fi

samtools sort -n -O BAM -o sorted.bam -T temp -@ 10 $bam
bedtools bamtofastq -i sorted.bam -fq $sample_1.fq -fq2 $sample_2.fq
rm sorted bam

bgzip $sample_1.fq
bgzip $sample_2.fq
