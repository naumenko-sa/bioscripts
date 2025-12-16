#!/bin/bash

# select variants with conflicing pathogenicity status
# install vembrane via conda and activate vembrane env
# if CLNSIG has a comma, it is parced as a tuple, i.e. "Pathogenic/Pathogenic,_low_penetrance|other" = ("Pathogenic/Pathogenic", "_low_penetrance|other")

# starts with Pathogenic or Likely_pathogetic, checksum matches
vembrane filter 'any(clnsig[:45] == "Conflicting_classifications_of_pathogenicity" for clnsig in INFO["CLNSIG"]) and (INFO["AF_EXAC"] > 0.001)' clinvar.vcf.gz > clinvar.conflicting.vcf

bgzip clinvar.conflicting.vcf
tabix clinvar.conflicting.vcf.gz

gunzip -c clinvar.conflicting.vcf.gz | grep "^#" > clinvar.conflicting.chr.vcf
gunzip -c clinvar.conflicting.vcf.gz | grep -v "^#" | awk '{print "chr"$0}' | sed s/chrMT/chrM/ >> clinvar.conflicting.chr.vcf
bgzip clinvar.conflicting.chr.vcf
tabix clinvar.conflicting.chr.vcf.gz

