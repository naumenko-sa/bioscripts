#!/bin/bash

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

#calculates read number in a bam file

samtools stats $bam | grep "^SN" | cut -f 2- | head -n1 | awk '{print $4}' > ${bam}.reads_number
