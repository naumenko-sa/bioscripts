#!/bin/bash

# $1 - input.hg19.vcf.gz
# $2 - output.vcf.gz
# $3 - hg19.fa
# $4 - hg38.fa
# $5 - chain.gz

bname=`basename $2 .vcf.gz`

# PRI tag has issues in Dragen
bcftools annotate -x FMT/PRI $1 | bcftools norm -m -both -f $3  | bgzip -c > $bname.norm.vcf.gz
tabix -p vcf $bname.vcf.gz

gatk LiftoverVcf --VERBOSITY ERROR \
  -I $bname.norm.vcf.gz \
  -O $bname.lifted.vcf.gz \
  -CHAIN $5 \
  -REJECT $bname.rejected.vcf.gz \
  -R $4

bcftools norm -f $4 -m -both $bname.lifted.vcf.gz | bgzip -c > $2
tabix -p vcf $2

rm $bname.lifted.vcf.gz $bname.lifted.vcf.gz.tbi
rm $bname.norm.vcf.gz $bname.norm.vcf.gz.tbi