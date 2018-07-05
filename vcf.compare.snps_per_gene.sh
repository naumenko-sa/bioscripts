#!/bin/bash

# writes table of variant hitting particular genes
# $1 = bed, $2 = vcf

#remove indels

bname=`basename $2 .vcf.gz`

bcftools filter -e 'TYPE="indel"' -Oz -o $bname.snps.vcf.gz  $2
tabix $bname.snps.vcf.gz

bedtools intersect -c -a $1 -b $bname.snps.vcf.gz > $2.gene_wise.table
for g in `cat $2.gene_wise.table | awk '{print $4}' | sort | uniq`
do 
    cat $2.gene_wise.table | grep $g | awk '{sum+=$5}END{print sum}' >> $2.variant_count
done
