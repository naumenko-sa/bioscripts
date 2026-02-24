#!/bin/bash

# select by AF > 0.001 after annotation - pathogenic variants that missed by freq arm
bname=`basename $1 .vcf.gz`
bcftools +fill-tags $1 -Oz -o $bname.max.vcf.gz -- -t 'max_af_calc:Float=fmax(max_af_regeneron|0, fmax(AFR_AF|0, AMR_AF|0))'
tabix $bname.max.vcf.gz
vembrane filter --annotation-key CSQ 'INFO["max_af_calc"] >= 0.001' $bname.max.vcf.gz > $bname.high_af_path.vcf
#  EAS_AF, EUR_AF, SAS_AF,
#                                                                gnomADe_AFR_AF, gnomADe_AMR_AF, gnomADe_EAS_AF,
#                                                                gnomADe_NFE_AF, gnomADe_SAS_AF, gnomADg_AFR_AF,
#                                                                gnomADg_AMR_AF, gnomADg_EAS_AF, gnomADg_MID_AF,
#                                                                gnomADg_NFE_AF, gnomADg_SAS_AF $1 > $bname.high_af_path.vcf

bgzip $bname.high_af_path.vcf
tabix $bname.high_af_path.vcf.gz


