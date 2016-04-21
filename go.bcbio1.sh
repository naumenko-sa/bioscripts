#!/bin/bash

#PBS -l walltime=50:00:00,nodes=1:ppn=16
#PBS -joe .
#PBS -d .
#PBS -l vmem=120g

echo "START: " `date`;

bcbio_nextgen.py ../config/project.yaml -n 16

echo "END: " `date`;


