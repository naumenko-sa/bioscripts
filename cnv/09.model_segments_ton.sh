#!/bin/bash

# https://gatk.broadinstitute.org/hc/en-us/articles/360042913531-ModelSegments
# run modelsegments without allelic counts - tumor only

# $1 = T.allelicCounts.tsv

bname=`basename $1 .allelicCounts.tsv`

gatk --java-options "-Xmx4g" ModelSegments \
--allelic-counts $1 \
--output . \
--output-prefix $bname
