#!/bin/bash

gatk VariantsToTable \
-V $1 \
-O $1.tsv \
-F CHROM -F POS -F ALT -F END -F SVLEN -F SVTYPE -F FOLD_CHANGE_LOG -F PROBES -F CN
