#!/bin/bash

# $1 = T.denoisedCR.tsv
# $2 hg38.dict

bname=`basename $1 .counts.denoisedCR.tsv`

gatk PlotModeledSegments \
--denoised-copy-ratios $1 \
--segments $bname.modelFinal.seg \
--sequence-dictionary $2 \
--minimum-contig-length 46709983 \
--output segment_plots \
--output-prefix $bname

