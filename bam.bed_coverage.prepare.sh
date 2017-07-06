#!/bin/bash

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=20g,mem=20g

sample=`echo $bam | sed s/.bam//`

samtools view -bh $bam $chr > $sample.$chr.bam
samtools index $sample.$chr.bam

