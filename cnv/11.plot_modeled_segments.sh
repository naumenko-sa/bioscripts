#!/bin/bash

# $1 = T.denoisedCR.tsv
# $2 hg38.dict

bname=`basename $1 .counts.denoisedCR.tsv`

gatk PlotModeledSegments \
--denoised-copy-ratios $1 \
--allelic-counts $bname.hets.tsv \
--segments $bname.modelFinal.seg \
--sequence-dictionary $2 \
--minimum-contig-length 10 \
--output segment_plots \
--output-prefix $bname

