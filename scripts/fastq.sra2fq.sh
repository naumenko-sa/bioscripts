#!/bin/bash

#PBS -l walltime=23:00:00,nodes=1:ppn=1
#PBS -d .
#PBS -joe
#PBS -l vmem=10g

[ -z $srr ] && srr=$1
[ -z $sample ] && sample=$2

module load sratoolkit
# download and convert in one go
# don't use -I - it causes bwa mem to failure as it requires identical read names

# also possible to fetch first - better on bad connections
# prefetch $srr, file goes to ~/ncbi/public/sra
fastq-dump.2.8.2 --gzip --split-files $srr
#or $sample.sra

mv ${srr}_1.fastq.gz ${sample}_1.fq.gz
mv ${srr}_2.fastq.gz ${sample}_2.fq.gz
