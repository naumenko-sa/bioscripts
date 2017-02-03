#!/bin/bash

module load gatk
java -Xmx4g -jar $GATK -R /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa \
     -T VariantFiltration -V $1 -o `echo $1 | sed s/vcf/DPLT10.vcf/` \
     --filterExpression " DP < 10" --filterName "DPLT10"