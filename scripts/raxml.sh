#!/bin/bash

#PBS -l walltime=23:00:00,nodes=1:ppn=5
#PBS -joe .
#PBS -d .
#PBS -l vmem=5g,mem=5g

if [ -z $align ]
then
    align=$1
fi

raxml=/hpf/largeprojects/ccmbio/naumenko/tools/standard-RAxML/raxmlHPC-PTHREADS-SSE3
$raxml -T 5 -s $align -n $align.tree -m GTRGAMMA --HKY85 -p 8593771075132634533
#-m PROTGAMMAWAG