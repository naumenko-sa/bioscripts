#!/bin/bash

# run Mutect2 in T only mode to produce calls for PON
# or TO calls without matched normals

# $1 = sample_N.bam
# #2 = panel.interval_list

# -tumor is deprecated

bname=`basename $1 .bam`

gatk Mutect2 \
-R /data/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa \
-I $1 \
-O $bname.vcf.gz \
--max-mnp-distance 0 \
--intervals $2 \
--interval-padding 50 \
--germline-resource af-only-gnomad.hg38.vcf.gz \
--genotype-germline-sites \
--native-pair-hmm-threads 16

tabix -f $bname.vcf.gz
