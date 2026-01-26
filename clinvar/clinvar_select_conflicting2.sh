#!/bin/bash

# select by AF > 0.001 after annotation - pathogenic variants that missed by freq arm
bname=`basename $1 .vcf.gz`
vembrane filter --annotation-key CSQ 'INFO["max_af_regeneron"] >= 0.001 or CSQ["AFR_AF"] >= 0.001 or CSQ["AMR_AF"] >= 0.001' > $bname.high_af_path.vcf
#  EAS_AF, EUR_AF, SAS_AF,
#                                                                gnomADe_AFR_AF, gnomADe_AMR_AF, gnomADe_EAS_AF,
#                                                                gnomADe_NFE_AF, gnomADe_SAS_AF, gnomADg_AFR_AF,
#                                                                gnomADg_AMR_AF, gnomADg_EAS_AF, gnomADg_MID_AF,
#                                                                gnomADg_NFE_AF, gnomADg_SAS_AF $1 > $bname.high_af_path.vcf

bgzip $bname.high_af_path.vcf
tabix $bname.high_af_path.vcf.gz


