#!/bin/bash

# select benign variants

# $1 clinvar.vcf.gz - nochr or clinvar.panel.vcf.gz - nochr

bname=`basename $1 .vcf.gz`
# add chr
gunzip -c $1 | grep "^#" | grep -v "contig=<ID=NT*" | sed -E 's/(contig=<ID=)([0-9]+)([^>]*>)/\1chr\2\3/g; s/(contig=<ID=)([XY])([^>]*>)/\1chr\2\3/g'  | sed s/"contig=<ID=MT"/"contig=<ID=chrM"/ > $bname.chr.vcf
gunzip -c $1 | grep -v "^#" | awk '{print "chr"$0}' | sed s/chrMT/chrM/ >> $bname.chr.vcf
bgzip $bname.chr.vcf
tabix $bname.chr.vcf.gz

gunzip -c $bname.chr.vcf.gz | grep "^#" > $bname.benign_conflicting.vcf
gunzip -c $bname.chr.vcf.gz | grep -v "^#" | grep "criteria_provided,_conflicting_classifications" | grep -E "(Benign|Likely_benign)" | grep -vE "CLNSIGCONF=(Pathogenic|Likely_pathogenic)" >> $bname.benign_conflicting.vcf
bgzip $bname.benign_conflicting.vcf
tabix $bname.benign_conflicting.vcf.gz

rm $bname.chr.vcf.*