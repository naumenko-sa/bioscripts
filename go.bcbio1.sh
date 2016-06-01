#!/bin/bash

#PBS -l walltime=500:00:00,nodes=1:ppn=20
#PBS -joe .
#PBS -d .
#PBS -l vmem=120g

echo "START: " `date`;

bcbio_nextgen.py ../config/${project}.yaml -n 20

echo "END: " `date`;


