#!/bin/bash

#SBATCH --job-name=bcbio
#SBATCH --mem=10G
#SBATCH --export=ALL
#SBATCH -t 7-50:00
#SBATCH -p core -n 1

# $1 = sort.bam
# $2 = interval.list - not gc annotated
# $3 = TSV for tsv output, default is hdf5

date

module load bcbio-nextgen/latest-devel

bname=`basename $1 .bam`

if [ -z $3 ]
then
    gatk --java-options '-Xms500m -Xmx131859m -XX:+UseSerialGC -Djava.io.tmpdir=.' \
    CollectReadCounts \
    -I $1 \
    -L $2 \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O $bname.counts.hdf5
else
   gatk --java-options '-Xms500m -Xmx131859m -XX:+UseSerialGC -Djava.io.tmpdir=.' \
    CollectReadCounts \
    -I $1 \
    -L $2 \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O $bname.counts.tsv \
    --format TSV
fi

date
