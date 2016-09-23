#!/bin/bash

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=50g,mem=50g

#module load gatk/3.6.0
#module load java/1.7.0

echo "START: " `date`;

module load java/1.8.0_91
#module load gatk
java -Xmx1g -version
which java

echo $project $threads

bcbio_nextgen.py ../config/${project}.yaml -n $threads

echo "END: " `date`;


