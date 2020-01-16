#!/bin/bash

gatk VariantsToTable \
-R ~/work/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa \
-V $1 \
-O $1.tsv \
-F CHROM -F POS -F REF -F ALT -F CSQ \
-U ALLOW_SEQ_DICT_INCOMPATIBILITY
