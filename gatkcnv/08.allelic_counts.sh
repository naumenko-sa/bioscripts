#!/bin/bash
# $1 = bam
# $2 = interval list
# $3 = reference

bname=`basename $1 .bam`

gatk --java-options "-Xmx3g" CollectAllelicCounts \
    -L $2 \
    -I $1 \
    -R $3 \
    -O $bname.allelicCounts.tsv
