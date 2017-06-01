#!/bin/bash

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=20g,mem=20g

module load java

cramtools bam -I $cram -O $bam -b -R $ref
