#!/bin/bash

#PBS -l walltime=230:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=20g,mem=20g

#lobSTR wrapper - calls short tandem duplications
#for single end reads

module load lobstr

PATH_TO_LOBSTR=/hpf/tools/centos6/lobstr/4.0.6/

lobSTR \
    --fastq \
    --gzip \
    -f ${sample}_1.fq.gz \
    --index-prefix $PATH_TO_LOBSTR/hg19_v3.0.2/lobstr_v3.0.2_hg19_ref/lobSTR_ -o $sample \
    --rg-sample $sample --rg-lib $sample

samtools sort $sample.aligned.bam -o $sample.sorted.bam -O BAM
samtools index $sample.sorted.bam

allelotype \
    --command classify \
    --bam $sample.sorted.bam \
    --noise_model $PATH_TO_LOBSTR/share/lobSTR/models/illumina_v3.pcrfree \
    --out $sample \
    --strinfo $PATH_TO_LOBSTR/hg19_v3.0.2/lobstr_v3.0.2_hg19_strinfo.tab \
    --index-prefix $PATH_TO_LOBSTR/hg19_v3.0.2/lobstr_v3.0.2_hg19_ref/lobSTR_
