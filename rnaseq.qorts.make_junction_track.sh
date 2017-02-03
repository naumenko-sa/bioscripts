#!/bin/bash

#PBS -l walltime=20:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=15g,mem=15g

#creates junction track for IGV browser from qorts output
#https://dl.dropboxusercontent.com/u/103621176/pipelineWalkthrough/example-walkthrough.pdf, page 27

#$1 = file = forJunctionSeq.txt.gz
#$2 = sample = output_name, sample.known.bed.gz

if [ -z $file ]
then
	file=$1
fi

if [ -z $sample ]
then
	sample=$2
fi

#1g is no enough
java -Xmx10g -jar ~/work/tools/bin/QoRTs.jar \
    makeJunctionTrack \
    --filenames $file \
    --trackTitle $sample \
    ~/work/tools/bcbio/genomes/Hsapiens/GRCh37/rnaseq/ref-transcripts.qorts.gff.gz \
    $sample.known.bed.gz
    