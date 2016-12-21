#!/bin/bash

#PBS -l walltime=1:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=15g,mem=15g

java -Xmx10g -jar ~/work/tools/bin/QoRTs.jar \
    mergeNovelSplices \
    --minCount 6 \
    outputData/countTables/ \
    outputData/sizeFactors.GEO.txt \
    inputData/annoFiles/anno.gtf.gz \
    outputData/countTables/