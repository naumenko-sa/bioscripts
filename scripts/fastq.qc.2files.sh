#!/bin/bash

if [ $# -lt "2" ]
then
    echo "Runs fastqc for left.fq [and right.fq]"
    echo "Usage: qc.sh numproc left.fq [right.fq]"
    exit 1
fi

fastqc -noextract -nogroup -t $1 $2 $3