#!/bin/bash

#PBS -l walltime=100:00:00,nodes=1:ppn=10
#PBS -joe .
#PBS -d .
#PBS -l vmem=120g

module load gatk/3.6.0

echo "START: " `date`;

bcbio_nextgen.py ../config/${project}.yaml -n 10

echo "END: " `date`;


