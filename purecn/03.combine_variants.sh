#!/bin/bash

vcf_files=""
for f in *.for_pon.vcf.gz
do 
    vcf_files="$vcf_files -V $f"
done

gatk3 -Xmx12g \
-T CombineVariants \
--minimumN 3 \
-R /data/genomes/Hsapiens/hg38/seq/hg38.fa \
-o snv_pon.vcf \
$vcf_files

bgzip snv_pon.vcf
tabix snv_pon.vcf.gz
