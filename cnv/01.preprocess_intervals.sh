#!/bin/bash

# $1 = coverage.bed

gatk PreprocessIntervals \
-R hg38.fa \
--interval-merging-rule OVERLAPPING_ONLY \
-O target.interval_list \
-L $1 \
--bin-length 0 \
--padding 250