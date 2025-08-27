#!/bin/bash

# convert from  GVCF to vcf and transfer filter values like PASS

# $1 - sample_name
# $2 - bed file
# $3 - reference

sample=$1
bed=$2
reference=$3

gatk GenotypeGVCFs \
      -R $reference \
      -V $sample.dragen.hard-filtered.gvcf.gz \
      -O $sample.dragen.hard-filtered.vcf.gz

tabix $sample.dragen.hard-filtered.vcf.gz

bedtools intersect -a $sample.dragen.hard-filtered.vcf.gz  \
                   -b $bed -header | \
         bgzip -c  > $sample.raw.vcf.gz

tabix $sample.raw.vcf.gz

bcftools annotate \
  -a $sample.dragen.hard-filtered.gvcf.gz \
  -c CHROM,POS,REF,FILTER \
  $sample.raw.vcf.gz \
  -O v  | bgzip -c  > $sample.vcf.gz
tabix $sample.vcf.gz

rm $sample.raw.vcf.gz $sample.raw.vcf.gz.tbi
rm $sample.dragen.hard-filtered.vcf.gz $sample.dragen.hard-filtered.vcf.gz.tbi
