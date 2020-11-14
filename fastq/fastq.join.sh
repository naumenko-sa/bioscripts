#!/bin/bash
#merges reads 

if [ $# -ne 4 ];
then
    echo "Usage : gofastqjoin.sh r1.fq r2.fq minimum-overlap output.fq";
    exit 1;
fi;

fastq-join $1 $2 -o $4 -m $3
wc -l ${4}join > ${4}join.wc