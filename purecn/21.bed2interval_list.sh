#!/bin/bash

# $1 = panel.bed
# $2 = hg38.fa
# out = panel.interval_list

bname=`basename $1 .bed`

gatk BedToIntervalList \
-I $1 \
-O $bname.interval_list \
-SD $bcbio/genomes/Hsapiens/hg38/seq/hg38.fa
#/data/bcbio/genomes/Hsapiens/hg38/seq/hg38.dict
