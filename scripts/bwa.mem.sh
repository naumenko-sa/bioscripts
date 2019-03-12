#!/bin/bash

#PBS -l walltime=20:00:00,nodes=1:ppn=10
#PBS -joe .
#PBS -d .
#PBS -l vmem=20g

bwa mem -t 10 SMN1.fasta $fq > $sample.SMN1.sam

samtools view -b $sample.SMN1.sam | samtools sort - > $sample.SMN1.bam
samtools index $sample.SMN1.bam
