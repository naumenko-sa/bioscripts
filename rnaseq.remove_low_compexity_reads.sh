#!/bin/bash

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g


if [ -z left ]
then
    left=$1
fi

if [ -z right ]
then
    right=$2
fi

if [ -z output ]
then
    output=$3
fi

gunzip $left $right
prinseq-lite.pl -fastq `echo $left | sed s/.gz//` -fastq2 `echo $right | sed s/.gz//` -log $output.prinseq.log -derep 12345 -lc_method dust -lc_threshold 7 -out_good $output