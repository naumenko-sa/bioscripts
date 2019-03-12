#!/bin/bash

if [ $# -lt 7 ];
then
    echo "Usage: gobowtie2.sh sample r1.fq r2.fq minins maxins threads orientation[rf|fr]";
    exit;
fi;

date
sample=$1 
f1=$2
f2=$3
output=${f1}_vs_${sample}
if [ ! -f $sample.build.log ];
then
    bowtie2-build $1 $sample > $sample.build.log;
fi;
bowtie2 -x $sample -1 $f1 -2 $f2  -S ${output}.sam --$7 --minins $4 --maxins $5 --threads $6 --no-unal --end-to-end &> ${output}.log
date

cat ${output}.sam | grep -v '^@' | awk '{print $2" "$9}' | grep -v 0$ | grep -v - | awk '{print $2}' > ${output}.is