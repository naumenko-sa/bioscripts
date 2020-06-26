#!/bin/bash

# $1 = T.counts.denoisedCR.tsv

bname=`basename $1 _T.counts.denoisedCR.tsv`

gatk --java-options "-Xmx4g" ModelSegments \
--denoised-copy-ratios $1 \
--allelic-counts ${bname}_T.allelicCounts.tsv \
--normal-allelic-counts ${bname}_N.allelicCounts.tsv \
--output . \
--output-prefix ${bname}_T
