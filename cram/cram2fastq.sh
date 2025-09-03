#!/bin/bash

# cram2bam
samtools view -O bam -o HG001_validation.bam HG001_validation.cram
# sort reads by name to output proper pairs
samtools sort -n -o HG001_validation.sorted.bam HG001_validation.bam
# bam2fastq
samtools fastq -1 HG001_validation_1.fq -2 HG001_validation_2.fq -s HG001_validation_singleton.fq HG001_validation.sorted.bam