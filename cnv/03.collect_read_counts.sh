#!/bin/bash

# $1 = sort.bam
# $2 = interval.list

# output TSV to plot graphs, or do hdf5

bname=`basename $1 -ready.bam`

gatk --java-options '-Xms500m -Xmx131859m -XX:+UseSerialGC -Djava.io.tmpdir=.' \
CollectReadCounts \
-I $1 \
-L $2 \
--interval-merging-rule OVERLAPPING_ONLY \
-O $bname.counts.tsv \
--format TSV
