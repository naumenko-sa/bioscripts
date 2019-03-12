#!/bin/bash

if [ $# -lt 3 ];
then
    echo "Usage: gobowtie2.sh sample r.fq threads";
    exit;
fi;

date
sample=$1 
f1=$2
output=${f1}_vs_${sample}
if [ ! -f $sample.build.log ];
then
    bowtie2-build $1 $sample > $sample.build.log;
fi;
bowtie2 -x $sample -U $f1 -S ${output}.sam --threads $3 --sensitive-local --no-unal &> ${output}.log
date

