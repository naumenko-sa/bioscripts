#!/bin/bash

#PBS -l walltime=1:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=15g,mem=15g

java -Xmx1g -jar ~/work/tools/bin/QoRTs.jar \
    makeJunctionTrack \
    --nonflatgtf \
    --filenames $1 \
    ../withNovel.forJunctionSeq.gff.gz \
    `echo $1 | sed s/txt.gz/known.bed.gz/`
