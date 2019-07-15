#!/bin/bash
# a stupid conversion of gvcf to vcf
# input gvcf is not gzipped
bname=`basename $1 .g.vcf`
cat $1 | grep "^#" > $bname.tmp.vcf
# also fixing FORMAT/AD which has extra 0 for NONREF
cat $1 | grep -v "^#" | awk '{if($5 != "<NON_REF>"){print $0}}' |  sed s/",<NON_REF>"// | sed s/",0:"/":"/ >> $bname.tmp.vcf
bgzip $bname.tmp.vcf
tabix $bname.tmp.vcf.gz
# remove FORMA/PL field
bcftools annotate -x FORMAT/PL $bname.tmp.vcf.gz -Oz > $bname.tmp2.vcf.gz
tabix $bname.tmp2.vcf.gz

~/cre/cre.vt.decompose.sh $bname.tmp2.vcf.gz

mv $bname.tmp2.decomposed.vcf.gz $bname.vcf.gz
tabix $bname.vcf.gz

rm $bname.tmp.* $bname.tmp2.*
