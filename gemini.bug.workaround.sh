#!/bin/bash

# there is a bug in gemini
# variant_impacts table is populated with the same values of hgvsc and hgvsp for a variant
# I have to extract this information from vcf file

if [ -z $family ]
then
    family=$1
fi

cat ${family}-ensemble.db.txt  | awk -F "\t" '{print $23"\t"$24"\t"$24}' | sed 1d | sed s/chr// > $family.bug.bed
bedtools intersect -header -a ${family}-ensemble.vcf.gz -b $family.bug.bed > $family.bug.vcf
gatk -R ~/work/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa -T VariantsToTable -V ${family}.bug.vcf -F CHROM -F POS -F REF -F ALT -F CSQ -o $family.bug.table
gemini.bug.parse_vep.pl $family.bug.table > $family.impacts.hgvs
