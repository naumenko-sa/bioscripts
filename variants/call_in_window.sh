#!/bin/bash

# call variants in a small window
# $1 - sample.bam
# $2 - window, chr:start-end

bname=`basename $1 .bam`

reference=~/work/grch37/GRCh37.fa.gz

samtools view -hb $1 $2 | samtools mpileup -E -uf $reference - | bcftools call -Ov -c -v - | bcftools view -i 'DP>10' > $bname.window.vcf
cat $bname.window.vcf  | awk -F '\t' '{print $1"\t"$2"\t"$4"\t"$5"\t"$NF}' > $bname.window.tsv
