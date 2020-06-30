#!/bin/bash

# $1 = sort.bam
# $2 = interval.list - not gc annotated
# $3 = TSV for tsv output, default is hdf5

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
