#!/bin/bash

#PBS -l walltime=10:00:00,nodes=1:ppn=10
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g

module load samtools
module load bedtools

samtools sort -n -O BAM -o sorted.bam -T temp -@ 10 $file
bedtools bamtofastq -i sorted.bam -fq ${pref}_1.fq -fq2 ${pref}_2.fq
