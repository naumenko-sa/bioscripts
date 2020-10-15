#!/bin/bash

# extract variants from mutect2 bcbio result vcf file
# select PASS variants
# filter on target variants

# $1 = bSample-mutect2-annotated.vcf.gz from bcbio

bname=`echo $1 | awk -F '-' '{print $1"."$2}' | awk '{print substr($0,2)}'`

echo "Sample:" $bname

bcftools view -f PASS -Oz $1 > $bname.pass.vcf.gz
tabix $bname.pass.vcf.gz

bedtools intersect \
-a $bname.pass.vcf.gz \
-b /n/shared_db/bcbio/biodata/genomes/Hsapiens/hg38/coverage/capture_regions/Exome-Agilent_V6.bed \
-header | bgzip > $bname.filtered.vcf.gz

tabix $bname.filtered.vcf.gz

gatk VariantsToTable \
-V $bname.filtered.vcf.gz \
-O $bname.tsv \
-F CHROM -F POS -F REF -F ALT -F ID -F DP -F ANN \
-GF GT -GF AD -GF AF

echo "Raw variants:" `gunzip -c $1 | grep -vc '^#'`
echo "PASS variants:" `gunzip -c $bname.pass.vcf.gz | grep -vc '^#'`
echo "On target variants:" `gunzip -c $bname.filtered.vcf.gz | grep -vc '^#'`

rm $bname.pass.vcf.gz
rm $bname.pass.vcf.gz.tbi
