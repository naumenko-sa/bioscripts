#!/bin/bash

# $1 = target.interval_list
# $2 = hg38.fa
# output = panel.gcannotated.tsv


bname=`basename $1 .interval_list`

gatk \
--java-options '-Xms500m -Xmx131859m -XX:+UseSerialGC -Djava.io.tmpdir=.' \
AnnotateIntervals \
-R $2 \
-L $1 \
--interval-merging-rule OVERLAPPING_ONLY \
-O $bname.gcannotated.tsv

