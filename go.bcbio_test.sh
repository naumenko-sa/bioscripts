#!/bin/bash

#PBS -l walltime=23:59:58,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=50g,mem=50g

echo "START: " `date`;

module load java/1.8.0_91
project=`cat families.txt | head -n $PBS_ARRAYID | tail -n1`

cd ${project}/work

#bcbio_nextgen.py ../config/${project}.yaml -n 5

echo "END: " `date`;


cd ../../