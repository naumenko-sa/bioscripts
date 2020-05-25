#!/bin/bash

# $1 = coverage.bed or could be interval_list file
# $2 = hg38.fa

# https://gatk.broadinstitute.org/hc/en-us/articles/360042914391-PreprocessIntervals
# recommended padding is 250 both sides - too much for exomes 100bp/exon? - 
# was optimal for TCGA exome data

bname=`basename $1 .bed`

gatk PreprocessIntervals \
-R $2 \
--interval-merging-rule OVERLAPPING_ONLY \
-O $bname.interval_list \
-L $1 \
--bin-length 0 \
--padding 0
