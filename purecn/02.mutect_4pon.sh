#!/bin/bash

# run Mutect2 in T only mode to produce calls for PON

bname=`basename $1 .bam`

gatk Mutect2 \
-R /data/genomes/Hsapiens/hg38/seq/hg38.fa \
-I $1 \
--max-mnp-distance 0 \
-O $bname.vcf.gz \
--intervals $2
