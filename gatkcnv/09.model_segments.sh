#!/bin/bash

# $1 = T.denoisedCR.tsv

bname=`echo $1 | awk -F "." '{print $1}'`

gatk --java-options "-Xmx4g" ModelSegments \
--denoised-copy-ratios $1 \
--output . \
--output-prefix $bname
