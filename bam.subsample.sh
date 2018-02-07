#!/bin/bash

#PBS -l walltime=10:00:00,nodes=1:ppn=10
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

if [ -z $bam ]
then
    bam=$1
fi

sambamba view -h -t 10 -s $fraction -f bam --subsampling-seed=123421 $bam -o `echo $bam | sed s/bam/subsample.bam/`
