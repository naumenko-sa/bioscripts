#!/bin/bash

# https://github.com/lindenb/jvarkit/wiki/BamStats04
# really fast: 2min for a big bam file and 100 genes
# arguments: bam and bed

#PBS -l walltime=2:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=21g

echo "bed: "$bed
echo "bam: "$bam

java -Xmx10G -jar ~/work/tools/jvarkit/dist/bamstats04.jar -B $bed $bam | sed s/"#"// | sed s/"\t"/","/g  > $bam.coverage.csv
