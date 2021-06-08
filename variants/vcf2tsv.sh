#!/bin/bash

gatk VariantsToTable \
-V $1 \
-O $1.tsv \
-F SAMPLE -F CHROM -F POS -F REF -F ALT -F DP -GF AF -F ANN
# CSQ AD
