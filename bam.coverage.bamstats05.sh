#!/bin/bash

# https://github.com/lindenb/jvarkit/wiki/BamStats05
# uses bamstats05 - groups coverage by gene
# really fast: 2min for a big bam file and 100 genes
# arguments: bam and bed

# protein coding genes are in ~/Desktop/reference_tables/protein_coding_genes.bed

#PBS -l walltime=2:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=21g

#wants java1.8
module load java

BAMSTATS_PATH=/hpf/largeprojects/ccmbio/naumenko/tools/jvarkit/dist
java -Xmx10G -jar ${BAMSTATS_PATH}/bamstats05.jar \
    -B=$bed $bam > $bam.coverage;
