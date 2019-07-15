#!/bin/bash

#PBS -l walltime=20:00:00,nodes=1:ppn=10
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

#samtools sort -n -O BAM -o sorted.bam -T temp -@ 10 $bam
bedtools bamtofastq -i $bam -fq ${sample}_1.fq -fq2 ${sample}_2.fq
#rm sorted.bam

bgzip ${sample}_1.fq
#bgzip ${sample}_2.fq
