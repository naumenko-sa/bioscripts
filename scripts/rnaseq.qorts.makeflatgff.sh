#!/bin/bash

#PBS -l walltime=1:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=15g,mem=15g

#module load java
java -Xmx10g -jar ~/work/tools/bin/QoRTs.jar makeFlatGff ref-transcripts.gtf ref-transcripts.qorts.gff
