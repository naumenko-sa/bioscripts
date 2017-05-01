#!/bin/bash

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g


#calculates features (reads) for RPKM calculation in R, outputs length of the genes

if [ -z $bam ]
then
    bam=$1
fi

featureCounts -T 8 \
	-g gene_id \
	-C \
	--largestOverlap \
	-a /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/rnaseq/ref-transcripts.gtf \
	-o $bam.rpkm_counts.txt $bam
