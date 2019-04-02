#!/bin/bash

#PBS -d .
#PBS -l walltime=10:00:00,nodes=1:ppn=24

#run with -o and without -o to get bootstrap support for all branches

if [ -z $align ]
then
    align=$1
fi

date
raxml=/hpf/largeprojects/ccmbio/naumenko/tools/standard-RAxML/raxmlHPC-PTHREADS-SSE3
$raxml -T 24 -s $align -n $align.tree -m GTRGAMMAI -p 8593771075132634533  -# 100
$raxml -T 24 -s $align -n $align.boot -m GTRGAMMAI -p 8593771075132634533 -b 1234567 -# 100
$raxml -T 24 -m GTRGAMMAI -p 8593771075132634533 -f b -t RAxML_bestTree.$align.tree -z RAxML_bootstrap.$align.boot -n $align.final
rm *tree.RUN*

#-f a - rapid bootstrap with the tree
#rapid bootstrapping
#-x 1234 
date
#-m PROTGAMMAWAG
