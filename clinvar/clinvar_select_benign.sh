#!/bin/bash

# select variants with conflicing pathogenicity status
# CLNSIG=criteria_provided,_conflicting_classifications
# no Pathogenic or Likely_pathogenic
# CLNSIGCONF=Uncertain_significance(7)|Benign(1)|Likely_benign(5)

# starts with Pathogenic or Likely_pathogetic, checksum matches
# AF_EXAC is not avaliable for some variants, so it is better to annotate with GNOMAD - don't consider AF_EXAC

# $1 clinvar.vcf.gz - nochr or clinvar.panel.vcf.gz - nochr
vembrane filter 'any(clnsig[:45] == "Conflicting_classifications_of_pathogenicity" for clnsig in INFO["CLNSIG"])' clinvar.vcf.gz > clinvar.conflicting.vcf

bgzip clinvar.conflicting.vcf
tabix clinvar.conflicting.vcf.gz

gunzip -c clinvar.conflicting.vcf.gz | grep "^#" | grep -v "contig=<ID=NT*" | sed -E 's/(contig=<ID=)([0-9]+)([^>]*>)/\1chr\2\3/g; s/(contig=<ID=)([XY])([^>]*>)/\1chr\2\3/g'  | sed s/"contig=<ID=MT"/"contig=<ID=chrM"/ > clinvar.conflicting.chr.vcf
gunzip -c clinvar.conflicting.vcf.gz | grep -v "^#" | awk '{print "chr"$0}' | sed s/chrMT/chrM/ >> clinvar.conflicting.chr.vcf
bgzip clinvar.conflicting.chr.vcf
tabix clinvar.conflicting.chr.vcf.gz

gunzip -c clinvar.conflicting.chr.vcf.gz | grep "^#" > clinvar.conflicting.pathogenic.vcf
gunzip -c clinvar.conflicting.chr.vcf.gz | egrep "CLNSIGCONF=(Pathogenic|Likely_pathogenic)" >> clinvar.conflicting.pathogenic.vcf
bgzip clinvar.conflicting.pathogenic.vcf
tabix clinvar.conflicting.pathogenic.vcf.gz
