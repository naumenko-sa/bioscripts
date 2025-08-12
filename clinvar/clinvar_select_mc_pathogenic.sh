#!/bin/bash

# select variants by MC (molecular consequence)
# get MC stat
# gunzip -c clinvar.PANELWES.missing.vcf.gz | grep -v "^#" |  awk -F 'MC=' '{print $2}' | awk -F ';' '{print $1}' | sort  | uniq -c | awk '{print $2"\t"$1}' | sort -k2,2nr > MC_stat.tsv

# don't report of MC is empty, usually these are deep upstream variants
# https://www.ncbi.nlm.nih.gov/snp/rs199469484

bname=`basename $1 .vcf.gz`

vembrane filter 'any(val in mc for val in ("frameshift_variant", "splice_donor_variant", "splice_acceptor_variant", \
                                           "missense_variant", "nonsense", "inframe_deletion", "initiator_codon_variant") for mc in INFO["MC"])' $1 > $bname.mc_coding.vcf
bgzip $bname.mc_coding.vcf
tabix $bname.mc_coding.vcf.gz

