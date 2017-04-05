#!/bin/bash

#https://github.com/lindenb/jvarkit/wiki/BamStats04
#uses bamstats05 - grouping coverage by gene
#arguments: bam and bed

#PBS -l walltime=2:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=21g

module load java
java -Xmx10G -jar ~/work/tools/jvarkit/dist/bamstats05.jar \
    -B=$bed $bam > $bam.coverage;
