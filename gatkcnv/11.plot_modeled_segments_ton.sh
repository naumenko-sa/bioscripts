#!/bin/bash

#SBATCH --job-name=bcbio
#SBATCH --mem=10G
#SBATCH --export=ALL
#SBATCH -t 7-50:00
#SBATCH -p core -n 1

# $1 = T.denoisedCR.tsv
# $2 hg38.dict

bname=`basename $1 .denoisedCR.tsv`

gatk PlotModeledSegments \
--denoised-copy-ratios $1 \
--segments $bname.modelFinal.seg \
--sequence-dictionary /projects/ngs/reference/UpdateGenomesBcbio/Hsapiens/hg38/seq/hg38.dict \
--minimum-contig-length 46709983 \
--output segment_plots \
--output-prefix $bname

