#!/bin/bash
# $1 = couts.standardizedCR.tsv
# $2 = hg38.dict

bname=`basename $1 .standardizedCR.tsv`

gatk PlotDenoisedCopyRatios \
--standardized-copy-ratios $1 \
--denoised-copy-ratios $bname.denoisedCR.tsv \
--sequence-dictionary /data/bcbio/genomes/Hsapiens/hg38/seq/hg38.dict \
--minimum-contig-length 46709983 \
--output plots \
--output-prefix $bname \
--point-size-copy-ratio 0.5
