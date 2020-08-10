#!/bin/bash

# $1 = panel.bed or could be interval_list file
# $2 = hg38.fa

# https://gatk.broadinstitute.org/hc/en-us/articles/360042914391-PreprocessIntervals
# recommended padding is 250 both sides - too much for exomes 100bp/exon? -
# was optimal for TCGA exome data

# some samples fail downstream when using 250 padding

# remove unplaced contigs and alts
# cat Exome-Agilent_V6.bed | grep -v random | grep -v alt | grep -v Un > panel.bed

. .profile

bname=`basename $1 .bed`

gatk PreprocessIntervals \
-R  $bcbio/genomes/Hsapiens/hg38/seq/hg38.fa \
--interval-merging-rule OVERLAPPING_ONLY \
-O $bname.padding250.interval_list \
-L $1 \
--bin-length 0 \
--padding 250
