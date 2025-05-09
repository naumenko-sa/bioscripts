#!/bin/bash

#$1 - bam or cram
#$2 - reference.fasta if cram

# has lots of messages, run as a job, or &> sample.out
# could be over 10G RAM for ~2G WES cram

echo "START:" `date`

file=$1
ext="${file#*.}"
sample="${file%%.*}"

echo $ext

if [[ "$ext" == "bam" ]];then
    gatk --java-options "-Xmx20G" MarkDuplicates -I $1 -O $sample.marked.bam -M $sample.dup_metrics.txt
    #--QUIET true --VERBOSITY ERROR
elif [[ "$ext" == "cram" ]];then
    # output is bam
    gatk --java-options "-Xmx20G" MarkDuplicates -I $1 -O $sample.marked.bam -M $sample.dup_metrics.txt -R $2
    #--QUIET true --VERBOSITY ERROR
fi

samtools index $sample.marked.bam

echo "END:" `date`
