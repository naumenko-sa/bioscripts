#!/bin/bash

# get variants missing in clinvar

# $1 - clinvar.pathogenic.vcf.gz or clinvar.selected.vcf.gz
# $2 - panel (should have panel.genes.chr.bed and panel.all_exons.chr.bed in the current dir)
# adds to the panel.summary.txt

panel=$2

# select variants which fall into gene boundaries
bedtools intersect -header -a $1 -b $panel.genes.bed | gzip -c > clinvar.$panel.vcf.gz

# select variants which are missing in the panel
bedtools intersect -header -v -a clinvar.$panel.vcf.gz -b $panel.all_exons.chr.bed | gzip -c > clinvar.$panel.missing.vcf.gz

echo "Missing pathogenic of Clinvar: " `gunzip -c clinvar.$panel.missing.vcf.gz | grep -vc "^#"` >> $panel.summary.txt
