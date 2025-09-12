#!/bin/bash

# $1 - calls
# $2 - truth
# $3 - bed

bname=`basename $1 .gvcf.gz`

# first get the BP resolution gVCF
gatk CombineGVCFs -R /data/reference/hg19_v2.fa \
--convert-to-base-pair-resolution true \
--call-genotypes true \
--variant $1 -O $bname.bp.gvcf.gz
tabix $bname.bp.gvcf.gz

# bedtools won't work with such vcf, unless fixed
# gunzip -c $bname.bp.gvcf.gz | grep "^#" > test.vcf
# gunzip -c $bname.bp.gvcf.gz | grep -v "^#" | awk -F '\t' '{if ($5 == "<NON_REF>") print $0}' | sed s/".\/."/"0\/0"/ | grep -E "(0/0)|(0|0)" | \
#     sed s/"GT:"/"GT\t:"/ | sed s/"0\/0:"/"0\/0\t:"/ | sed s/"0|0"/"0|0\t:"/ | cut -f 1-9,11 | sed s/"<NON_REF>"/"."/ >> test.vcf
# bgzip test.vcf
# tabix test.vcf.gz


rtg vcfeval \
  -b $2 \
  -c $bname.bp.gvcf.gz \
  -t /data/reference/hg19_v2.SDF \
  -o eval_annotated \
  --bed-regions $3 \
  --output-mode annotate \
  --all-records
 
echo "TN: " `gunzip -c eval_annotated/calls.vcf.gz | grep -c "CALL=IGN"`

# no +1 because bed is [0-based, 1-based) [10,20) > [11,20] or [11, 21)
echo "Total BP in bed:" `cat $3 | awk -F '\t' '{sum+=$3-$2}END{print sum}'`
