#!/bin/bash

#PBS -l walltime=100:00:00,nodes=1:ppn=12
#PBS -d .

if [ $# -lt "2" ]
then
    echo "Runs fastqc for left.fq [and right.fq]"
    echo "Usage: qc.sh numproc left.fq [right.fq]"
    exit 1
fi

fastqc -noextract -nogroup -t $1 $2 $3