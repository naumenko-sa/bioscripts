#!/bin/bash
# also vt rminfo but it requires exact case matching
bname=`basename $1 .vcf.gz`
bcftools annotate -x ^INFO/AC,INFO/AF,INFO/AN,INFO/BaseQRankSum,INFO/ClippingRankSum,INFO/DB,INFO/DP,INFO/ExcessHet,INFO/FS,INFO/MLEAC,INFO/MLEAF,INFO/MQ,INFO/MQ0,INFO/MQRankSum,INFO/QD,INFO/ReadPosRankSum,INFO/SOR -Oz -o $bname.no_anno.vcf.gz $1
tabix $bname.no_anno.vcf.gz
