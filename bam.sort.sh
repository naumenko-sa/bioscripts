#!/bin/bash

# faster to sort on computing nodes

#PBS -d .

if [ -z $bam ]
then
    bam=$1
fi

# sam2bam
#samtools view -bS $1 > $bamname

sortedbam=`echo $1 | sed s/bam/sorted.bam/`
samtools sort $bam $sortedbam
samtools index $sortedbam

