#!/bin/bash

# select pathogenic variants reviewed by expert panel or supported by multiple labs
# install vembrane via conda and activate vembrane env
# if CLNSIG has a comma, it is parced as a tuple, i.e. "Pathogenic/Pathogenic,_low_penetrance|other" = ("Pathogenic/Pathogenic", "_low_penetrance|other")

# starts with Pathogenic or Likely_pathogetic, checksum matches
vembrane filter 'any(clnsig[:10]=="Pathogenic" or clnsig[:17]=="Likely_pathogenic" for clnsig in INFO["CLNSIG"])' clinvar.vcf.gz > clinvar.pathogenic.vcf

# this does not work for the bottom 4 and P/LP/P,_low_penetrance
# vembrane filter 'any(clnsig in INFO["CLNSIG"] for clnsig in ("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic", \
#                     "Pathogenic|other", "Pathogenic|drug_response", "Likely_pathogenic|drug_response", "Pathogenic|risk_factor", \
#                     "Likely_pathogenic/Likely_risk_allele", "Pathogenic/Likely_risk_allele", "Likely_pathogenic|risk_factor", \
#                     "Pathogenic|Affects", "Likely_pathogenic,_low_penetrance", "Pathogenic/Likely_pathogenic|other", "Likely_pathogenic|association", \
#                     "Pathogenic/Likely_pathogenic/Likely_risk_allele", \
#                     "Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance", \
#                     "Pathogenic|association", "Pathogenic/Likely_pathogenic|risk_factor", "Likely_pathogenic|Affects", "Pathogenic|protective", \
#                     "Likely_pathogenic|other", "Pathogenic|association|protective", "Pathogenic|confers_sensitivity",  \
#                     "Pathogenic/Likely_pathogenic|association", \
#                     "Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance|other", \
#                     "Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance|risk_factor", \
#                     "Pathogenic/Pathogenic,_low_penetrance|other", \
#                     "Pathogenic/Pathogenic,_low_penetrance|other|risk_factor"))' clinvar.vcf.gz > clinvar.pathogenic.vcf

bgzip clinvar.pathogenic.vcf
tabix clinvar.pathogenic.vcf.gz

gunzip -c clinvar.pathogenic.vcf.gz | grep "^#" > clinvar.pathogenic.chr.vcf
gunzip -c clinvar.pathogenic.vcf.gz | grep -v "^#" | awk '{print "chr"$0}' | sed s/chrMT/chrM/ >> clinvar.pathogenic.chr.vcf
bgzip clinvar.pathonetic.chr.vcf
tabix clinvar.pathogenic.chr.vcf.gz

# same story with CLNREVSTAT, when there is a comma, it parses as a tuple
# we need to extract
# criteria_provided,_multiple_submitters,_no_conflicts
# reviewed_by_expert_panel

vembrane filter '(INFO["CLNREVSTAT"][0] == "criteria_provided" and INFO["CLNREVSTAT"][1] == "_multiple_submitters" and INFO["CLNREVSTAT"][2] == "_no_conflicts") or \
                 (INFO["CLNREVSTAT"][0] == "reviewed_by_expert_panel")' clinvar.pathogenic.vcf.gz > clinvar.selected.vcf

bgzip clinvar.selected.vcf
tabix clinvar.selected.vcf.gz

# clinvar is nochr

gunzip -c clinvar.selected.vcf.gz | grep "^#" > clinvar.selected.chr.vcf
gunzip -c clinvar.selected.vcf.gz | grep -v "^#" | awk '{print "chr"$0}' | sed s/chrMT/chrM/ >> clinvar.selected.chr.vcf
bgzip clinvar.selected.chr.vcf
tabix clinvar.selected.chr.vcf.gz
