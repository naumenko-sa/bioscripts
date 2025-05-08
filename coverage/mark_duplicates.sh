#!/bin/bash

sample=`basename $1 .bam`

gatk MarkDuplicates -I $1 -O $sample.marked.bam -M $sample.dup_metrics.txt
samtools index $sample.marked.bam
