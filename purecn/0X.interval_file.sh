#!/bin/bash

# $1 = panel.bed
# $2 = GCA_000001405.15_GRCh38_no_alt_analysis_set_100.bw

PURECN=/data/bcbio/anaconda/envs/r36/lib/R/library/PureCN/extdata

export PATH=/data/bcbio/anaconda/envs/r36/bin:$PATH

which Rscript

Rscript $PURECN/IntervalFile.R \
--infile $1 \
--fasta /data/genomes/Hsapiens/hg38/seq/hg38.fa \
--outfile intvervals.txt \
--offtarget \
--genome hg38 \
--export baits_optimized_hg38.bed \
--mappability $2
