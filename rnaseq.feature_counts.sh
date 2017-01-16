#!/bin/bash

#PBS -l walltime=2:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

if [ -z bam ];
then
    bam=$1
fi
    

featureCounts  \
	-T 8 \ #threads
	-F GTF \
	-t exon \
	-g gene_id \ 
	-s 0 
	-C \
	--largestOverlap \
	-a /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/rnaseq/ref-transcripts.gtf \
	-o $bam.counts.txt \
	$bam
