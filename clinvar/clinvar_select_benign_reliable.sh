#!/bin/bash

# $1 - clinvar.vcf.gz or clinvar.panel.vcf.gz, no chr

# select benign  variants reviewed by expert panel or supported by multiple labs
# if CLNSIG has a comma, it is parced as a tuple, i.e. "Pathogenic/Pathogenic,_low_penetrance|other" = ("Pathogenic/Pathogenic", "_low_penetrance|other")

bname=`basename $1 .vcf.gz`

# add chr
gunzip -c $1 | grep "^#" | grep -v "contig=<ID=NT*" | sed -E 's/(contig=<ID=)([0-9]+)([^>]*>)/\1chr\2\3/g; s/(contig=<ID=)([XY])([^>]*>)/\1chr\2\3/g'  | sed s/"contig=<ID=MT"/"contig=<ID=chrM"/ > $bname.chr.vcf
gunzip -c $1 | grep -v "^#" | awk '{print "chr"$0}' | sed s/chrMT/chrM/ >> $bname.chr.vcf
bgzip $bname.chr.vcf
tabix $bname.chr.vcf.gz

# select benign and likely benign
# this is done, first otherwise vembrane errors out
vembrane filter 'any(clnsig[:6]=="Benign" or clnsig[:13]=="Likely_benign" for clnsig in INFO["CLNSIG"])' $bname.chr.vcf.gz > $bname.benign.vcf
bgzip $bname.benign.vcf
tabix $bname.benign.vcf.gz

vembrane filter '(INFO["CLNREVSTAT"][0] == "criteria_provided" and INFO["CLNREVSTAT"][1] == "_multiple_submitters" and INFO["CLNREVSTAT"][2] == "_no_conflicts") or \
                 (INFO["CLNREVSTAT"][0] == "reviewed_by_expert_panel") or \
                 (INFO["CLNREVSTAT"][0] == "practice_guideline")' $bname.benign.vcf.gz > $bname.benign_reliable.vcf
bgzip $bname.benign_reliable.vcf
tabix $bname.benign_reliable.vcf.gz

rm $bname.chr.vcf* $bname.benign.*