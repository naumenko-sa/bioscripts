#!/bin/bash

# $1 - input.vcf (not gz)
# $2 - hg19.fa
# $3 - hg38.fa
# $4 - chain.gz

bname=`basename $1 .vcf`

# the SNV liftover is not working for CNV
# convert to bed
# liftoverbed
# reconstruct VCF

