#!/bin/bash

# for clinvar
# activate vembrane conda environment

bname=`basename $1 .vcf.gz`

vembrane table 'INFO["RS"], CHROM, POS, REF, ALT, INFO["AF_EXAC"], INFO["GENEINFO"], INFO["MC"], INFO["CLNSIGCONF"], INFO["ALLELEID"]' $1 > $bname.tsv
