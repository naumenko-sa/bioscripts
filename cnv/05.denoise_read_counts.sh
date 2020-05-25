#!/bin/bash
# $1 = input bam.counts.tsv
# $2 = pon.hdf5

bname=`basename $1 .tsv`

gatk --java-options "-Xmx12g" \
DenoiseReadCounts \
-I $1 \
--count-panel-of-normals $2 \
--standardized-copy-ratios $bname.standardizedCR.tsv \
--denoised-copy-ratios $bname.denoisedCR.tsv
