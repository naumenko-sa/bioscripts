#!/bin/bash

# $1 - input.vcf

bname=`basename $1 .vcf`

bcftools view -H $1 | \
  awk -F'\t' '{
    chrom=$1; pos=$2; id=$3; ref=$4; alt=$5; info=$8;
    end = pos;
    if (match(info, /END=([0-9]+)/, a)) end=a[1];
    start = pos-1; print chrom"\t"start"\t"end"\t"id"\t0\t+\t"ref"\t"alt"\t"info
  }' > $bname.bed
