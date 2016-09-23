#!/bin/bash

#convert vep.vcf annotations for export to xls


#cat $1 | grep -v "^#" | grep PASS | awk '{for(i=1;i<8;i++)printf $i"\t";printf $9"\t"$10;gsub(";CSQ","\nCSQ",$8);gsub(",","\nCSQ=",$8);print $8}' | sed s/"|"/"\t"/g > $1.txt
#"Consequence annotations from Ensembl VEP. Format: Consequence|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|CANONICAL|CCDS|LoF|LoF_filter|LoF_flags">

echo -e "CHR\tPOS\tID\tREF\tALT\tConsequence\tCodons\tAmino_acids\tGene\tSYMBOL\tFeature\tEXON\tPolyPhen\tSIFT\tProtein_position\tBIOTYPE\tCANONICAL\tCCDS\tLoF\tLoF_filter\tLoF_flags" > $1.txt

cat $1 | grep -v "^#"  | grep PASS |  sed s/";CSQ="/"\t"/ | cut -f 1-5,9 | awk -F ';' '{print $1}' | awk '{n=split($6,t,",");for(i=1;i<n;++i){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"t[i]}}' | sed s/"|"/"\t"/g >> $1.txt
