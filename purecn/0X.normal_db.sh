#!/bin/bash

# $1 = snv_pon.vcf.gz from Mutect2 PON
PURECN=/data/bcbio/anaconda/envs/r36/lib/R/library/PureCN/extdata

export PATH=/data/bcbio/anaconda/envs/r36/bin:$PATH

which Rscript

Rscript $PURECN/NormalDB.R --help

#Rscript $PURECN/NormalDB.R \
#--outdir . \
#--normal_panel $1 \
#--assay exome_idt_v1 \
#--genome hg38 \
#--force
--coveragefiles
