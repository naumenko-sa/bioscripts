#!/bin/bash

# $1 = coverage.bed
# $2 = hg38.fa

bname=`basename $1 .bed`

gatk PreprocessIntervals \
-R $2 \
--interval-merging-rule OVERLAPPING_ONLY \
-O $bname.interval_list \
-L $1 \
--bin-length 0 \
--padding 250
