#!/bin/bash

#PBS -l walltime=23:00:00,nodes=1:ppn=1
#PBS -d .
#PBS -joe
#PBS -l vmem=10g

[ -z $srr ] && srr=$1
[ -z $sample ] && sample=$2

#2.8.2 on data, 2.8.0 on compute
module load sratoolkit
# download and convert in one go
# don't use -I - it causes bwa mem to failure as it requires identical read names

# also possible to fetch first - better on bad connections
# prefetch $srr, file goes to ~/ncbi/public/sra
# note srr is SRRNNNN without .sra
fastq-dump --gzip --split-files $srr.sra
#or $sample.sra

mv ${srr}_1.fastq.gz ${sample}_1.fq.gz
mv ${srr}_2.fastq.gz ${sample}_2.fq.gz
