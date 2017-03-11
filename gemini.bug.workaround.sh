#!/bin/bash

# there is a bug in gemini
# variant_impacts table is populated with the same values of hgvsc and hgvsp for a variant
# I have to extract this information from vcf file

if [ -z $family ]
then
    family=$1
fi

#extract positions of the variants in the final report
cat ${family}.no_vep.decomposed.vepeffects.db.txt  | awk -F "\t" '{print $23"\t"$24}' | sed 1d | sed s/chr// > $family.bug.tab
#extract variants with their VEP annotations
bcftools view -o $family.bug.tosort.vcf -R $family.bug.tab ${family}.no_vep.decomposed.vepeffects_refseq.vcf.gz 
#sort vcf
picard SortVcf I=$family.bug.tosort.vcf O=$family.bug.vcf
#convert vcf to table
gatk -R ~/work/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa -T VariantsToTable -V ${family}.bug.vcf -F CHROM -F POS -F REF -F ALT -F CSQ -o $family.bug.table -U ALLOW_SEQ_DICT_INCOMPATIBILITY
#parse vcf and VEP CSQ field
gemini.bug.parse_vep.pl $family.bug.table > $family.ref_seq.hgvs
