#!/bin/bash

# $1 = panel.interval_list
# output = panel.gcannotated.tsv

. .profile

bname=`basename $1 .interval_list`

gatk \
--java-options '-Xms500m -Xmx131859m -XX:+UseSerialGC -Djava.io.tmpdir=.' \
AnnotateIntervals \
-R  $bcbio/genomes/Hsapiens/hg38/seq/hg38.fa \
-L $1 \
--interval-merging-rule OVERLAPPING_ONLY \
-O $bname.gcannotated.tsv
