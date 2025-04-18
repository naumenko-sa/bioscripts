#!/bin/bash

# maybe because NUMBER=. for CLNREVSTAT it is a tuple not string

vembrane filter 'any(clnsig in INFO["CLNSIG"] for clnsig in ("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic", \
                                                             "Pathogenic|drug_response", "Likely_pathogenic|drug_response")) and \
                 ("reviewed_by_expert_panel" in INFO["CLNREVSTAT"])' clinvar.vcf.gz > clinvar.selected.vcf
bgzip clinvar.selected.vcf
tabix clinvar.selected.vcf.gz

# clinvar is nochr

gunzip -c clinvar.selected.vcf.gz | grep "^#" > clinvar.selected.chr.vcf
gunzip -c clinvar.selected.vcf.gz | grep -v "^#" | awk '{print "chr"$0}' | sed s/chrMT/chrM/ >> clinvar.selected.chr.vcf
bgzip clinvar.selected.chr.vcf
tabix clinvar.selected.chr.vcf.gz

