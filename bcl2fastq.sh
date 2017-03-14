#!/bin/bash

#PBS -l walltime=23:00:00,nodes=1:ppn=10
#PBS -joe .
#PBS -d .
#PBS -l vmem=20g,mem=20g

if [ -z $sample_sheet ];
then
    sample_sheet=$1
fi

#2.19 needs 20G of RAM for 10 threads
module load bcl2fastq

#I use no-lane splitting, and the result fastq files could not be uploaded to basespace
bcl2fastq --sample-sheet $sample_sheet -p  10 --no-lane-splitting
