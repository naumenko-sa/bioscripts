#!/bin/bash

#PBS -l walltime=359:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=50g,mem=50g

echo "START: " `date`;

module load java/1.8.0_91
echo $project $threads

bcbio_nextgen.py ../config/${project}.yaml -n $threads

echo "END: " `date`;


