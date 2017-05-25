#!/bin/bash

#PBS -l walltime=10:00:00,nodes=1:ppn=10
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g

if [ -z $bam ]
then
    bam=$1
fi

if [ -z $pref ]
then
    pref=$2
fi

samtools sort -n -O BAM -o sorted.bam -T temp -@ 10 $bam
bedtools bamtofastq -i sorted.bam -fq ${pref}_1.fq -fq2 ${pref}_2.fq

bgzip ${pref}_1.fq
bgzip ${pref}_2.fq
