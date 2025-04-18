#!/bin/bash

gunzip -c clinvar.pathogenic.sorted.vcf.gz | grep "^#" > clinvar.pathogenic.panela.vcf
bedtools intersect -a clinvar.pathogenic.sorted.vcf.gz -b ../coverage/panela/panela.genes.merged.bed >> clinvar.pathogenic.panela.vcf
bgzip clinvar.pathogenic.panela.vcf 
tabix clinvar.pathogenic.panela.vcf.gz 