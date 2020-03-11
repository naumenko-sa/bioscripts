#!/bin/bash

# bwa alignment of collapsed umi with alt

#SBATCH --job-name=bcbio
#SBATCH --mem=96G
#SBATCH --export=ALL
#SBATCH -t 2-00:00
#SBATCH -p core -n 48

bcbio_path=

bwa mem  -C -c 250 -M -t 48 \
-R '@RG\tID:062_01_E01_14_3\tPL:illumina\tPU:062_01_E01_14_3\tSM:062_01_E01_14_3' \
-v 1 /projects/ngs/reference/genomes/Hsapiens/hg38/bwa/hg38.fa \
062-01_E01_14_3-cumi-1.fq.gz \
062-01_E01_14_3-cumi-2.fq.gz | \
k8 ${bcbio_path/anaconda/share/bwakit-0.7.15-1/bwa-postalt.js \
-p 062_01_E01_14_3-sort-cumi.bam.hla /projects/ngs/reference/genomes/Hsapiens/hg38/bwa/hg38.fa.alt | \
samtools sort -@ 48 -m 2G  -T sort-cumi-sorttmp -o 062-01_E01_14_3-sort-cumi.bam /dev/stdin
