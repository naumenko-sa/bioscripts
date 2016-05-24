#!/bin/bash

#PBS -l walltime=100:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=20g

echo "START: " `date`;

bcbio_nextgen.py ../config/${project}.yaml -n 1

echo "END: " `date`;


