#!/bin/bash

# select variants by MC (molecular consequence)
# get MC stat
# gunzip -c clinvar.PANELWES.missing.vcf.gz | grep -v "^#" |  awk -F 'MC=' '{print $2}' | awk -F ';' '{print $1}' | sort  | uniq -c | awk '{print $2"\t"$1}' | sort -k2,2nr > MC_stat.tsv

# don't report of MC is empty, usually these are deep upstream variants
# https://www.ncbi.nlm.nih.gov/snp/rs199469484

bname=`basename $1 .vcf.gz`

vembrane filter 'any("frameshift_variant" in mc or \
                     "splice_donor_variant" in mc or \
                     "missense_variant" in mc or \
                     "nonsense" in mc or \
                     "splice_acceptor_variant" in mc or \
                     "inframe_deletion" in mc or \
                     "initiator_codon_variant" in mc \
                     for mc in INFO["MC"])' > $bname.mc_coding.vcf
bgzip $bname.mc_coding.vcf
tabix $bname.mc_coding.vcf.gz

