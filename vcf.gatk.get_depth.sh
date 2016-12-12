#!/bin/bash

bname=`echo $1 | sed s/.vcf.gz//`

gatk -R ~/work/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa -T VariantsToTable -V $1 \
     -F CHROM -F POS -F REF -F ALT -F DP -GF DP -GF AD -o $bname.table --allowMissingData 
