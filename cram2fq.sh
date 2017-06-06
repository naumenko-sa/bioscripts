#!/bin/bash

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=20g,mem=20g

# reference should be decompressed and indexed: using one from bcbio
# usually crams come with hg19 reference (chr), not GRCh37


module load java

cramtools fastq -Xmx10g -F $sample --skip-md5-check \
		-z \
		-I $cram \
		-R /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/hg19/seq/hg19.fa
