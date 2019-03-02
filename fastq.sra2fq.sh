#!/bin/bash

#PBS -l walltime=23:00:00,nodes=1:ppn=1
#PBS -d .
#PBS -joe
#PBS -l vmem=10g

[ -z $srr ] && srr=$1
[ -z $sample ] && sample=$2

module load sratoolkit
#download and convert in one go
fastq-dump -I --gzip --split-files $srr
#$sample.sra

mv ${srr}_1.fastq.gz ${sample}_1.fq.gz
mv ${srr}_2.fastq.gz ${sample}_2.fq.gz
