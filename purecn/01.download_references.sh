#!/bin/bash -l

. .profile
wget -c https://s3.amazonaws.com/purecn/GCA_000001405.15_GRCh38_no_alt_analysis_set_100.bw

# hg38_simpleRepeats.bed
Rscript 01.UCSC_repeats.R

# copy 2 files to /bcbio/anaconda/envs/r36/lib/R/library/PureCN/extdata