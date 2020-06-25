#!/bin/bash

# run Mutect2 in T only mode to produce calls for PON
# $1 = sample_N.bam
# #2 = panel.interval_list

bname=`basename $1 .bam`

gatk Mutect2 \
-R /bcbio/genomes/Hsapiens/hg38/seq/hg38.fa \
-I $1 \
-O $bname.for_pon.vcf.gz \
-tumor S1_N \
--max-mnp-distance 0 \
--intervals $2 \
--interval-padding 50 \
--germline-resoure af-only-gnomad.hg38.vcf.gz

tabix $bname.for_pon.vcf.gz
