#!/bin/bash

#https://github.com/lindenb/jvarkit/wiki/BamStats04

#PBS -l walltime=2:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=21g

module load java
java -Xmx10G -jar ~/work/tools/jvarkit/dist/bamstats04.jar \
    BED=$bed I=$bam > $bam.coverage;
