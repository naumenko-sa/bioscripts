#!/bin/bash

clinvar_select_benign_reliable.sh clinvar.panel.vcf.gz
clinvar_select_benign_conflicting.sh clinvar.panel.vcf.gz
bcftools merge -m none clinvar.panel.benign_reliable.vcf.gz clinvar.panel.benign_conflicting.vcf.gz | bcftools sort -Oz -o clinvar.panel.benign.vcf.gz
tabix clinvar.panel.benign.vcf.gz

gunzip -c clinvar.panel.benign.vcf.gz | grep -v "^#" | cut -f 1,2,4,5 > clinvar.benign.tsv
