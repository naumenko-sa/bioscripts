#!/bin/bash

# $1 = bam
# $2 = interval list

bname=`basename $1 -ready.bam`

gatk --java-options "-Xmx3g" CollectAllelicCounts \
    -L $2 \
    -I $1 \
    -R $3 \
    -O $bname.allelicCounts.tsv
