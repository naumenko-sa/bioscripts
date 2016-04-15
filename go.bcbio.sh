#!/bin/bash

#PBS -l walltime=50:00:00,nodes=1:ppn=16
#PBS -joe .
#PBS -d .
#PBS -l vmem=80g

# exceeded 47 GB on test exome data

bcbio_nextgen.py ../config/NA12878-exome-methodcmp.yaml -n 16

echo "END: " `date`;


