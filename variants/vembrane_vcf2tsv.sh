#!/bin/bash

# activate vembrane conda environment

bname=`basename $1 .vcf.gz`

vembrane table --annotation-key CSQ 'INFO["RS"], CHROM, POS, REF, ALT, INFO["ALLELEID"], INFO["CLNDN"], INFO["CLNHGVS"], INFO["GENEINFO"], INFO["CLNVI"]' $1 > $bname.tsv
