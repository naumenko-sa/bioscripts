#!/bin/bash

PURECN=/data/bcbio/anaconda/envs/r36/lib/R/library/PureCN/extdata

export PATH=/data/bcbio/anaconda/envs/r36/bin:$PATH

#which Rscript

# $1 = sample_id
# $2 = mapping_bias.rds
Rscript $PURECN/PureCN.R \
--sampleid $1 \
--out $1.out \
--tumor ${1}-target-coverage.hdf5 \
--logratiofile ${1}-crdenoised.tsv \
--segfile $1.modelFinal.seg \
--mappingbiasfile $2 \
--vcf $1.vcf.gz \
--statsfile $1.vcf.gz.stats \
--genome hg38 \
--funsegmentation Hclust \
--force \
--postoptimize \
--seed 123
