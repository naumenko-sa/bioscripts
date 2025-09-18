#!/bin/bash

# validate a vcf file against genome in a bottle calls for NA12878
# $1 - file.vcf.gz
# $2 - truth.vcf.gz
# $3 - regions.bed
# $4 - reference.SDF

vcf=$1
truth=$2
bed=$3
#reference.SDF
reference=$4


bname=`basename $1 .vcf.gz`

# rtg manual
# https://github.com/RealTimeGenomics/rtg-tools/blob/master/installer/resources/tools/RTGOperationsManual.pdf

# bedtools intersect -nonamecheck -a NA12878-sort-callable_sample.bed -b GiaB_v2_19_regions.bed > NA12878-sort-callable_sample-NA12878-wrm.bed

# uses PASS variants only
export RTG_JAVA_OPTS='-Xms750m' && export RTG_MEM=9100m && \
   rtg vcfeval --threads 5 -b $truth \
   --bed-regions $bed \
   -c $vcf \
   -t $reference \
   -o ${bname}_rtg --vcf-score-field='GQ' 
#   --all-records

for f in {tp-baseline,fp,fn}
do
    echo snp $f `bcftools view --types snps ${bname}_rtg/${f}.vcf.gz | grep -vc "^#"` >> $1.stat
    echo indels $f `bcftools view --exclude-types snps ${bname}_rtg/${f}.vcf.gz | grep -vc "^#"` >> $1.stat
done

