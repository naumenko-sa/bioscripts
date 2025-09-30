#!/bin/bash

# $1 - input.vcf (not gz)
# $2 - hg19.fa
# $3 - hg38.fa
# $4 - chain.gz

bname=`basename $1 .vcf`

bcftools norm -m -both -f $2 $1 | bgzip -c > $bname.norm.vcf.gz
tabix -p vcf $bname.vcf.gz

gatk LiftoverVcf \
  I=$bname.norm.vcf.gz \
  O=$bname.lifted.vcf.gz \
  CHAIN=$4 \
  REJECT=$bname.rejected.vcf.gz \
  R=$3

bcftools norm -f $3 -m -both $bname.lifted.vcf.gz | bgzip -c > $bname.lifted.norm.vcf.gz
tabix -p vcf $bname.lifted.norm.vcf.gz
