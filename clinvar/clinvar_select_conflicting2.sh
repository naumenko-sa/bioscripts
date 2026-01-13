#!/bin/bash

# select by AF > 0.001 after annotation - pathogenic variants that missed by freq arm
bname=`basename $1 .vcf.gz`
vembrane filter 'INFO["max_af_regeneron"] >= 0.001 or INFO["max_af_regeneron"] is NA' $1 > $bname.high_af_path.vcf

bgzip $bname.high_af_path.vcf
tabix $bname.high_af_path.vcf.gz


