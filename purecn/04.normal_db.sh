#!/bin/bash


PURECN=/data/bcbio/anaconda/envs/r36/lib/R/library/PureCN/extdata

export PATH=/data/bcbio/anaconda/envs/r36/bin:$PATH

which Rscript

Rscript $PURECN/NormalDB.R --help

#Rscript $PURECN/NormalDB.R \
#--outdir . \
#--normal_panel \
#--assay assayn \
#--genome hg38 \
#--force
