#/bin/bash

while read chr pos;
do 
    tabix mgp.v5.merged.snps_all.dbSNP142.vcf.gz $chr:$pos-$pos >> $1.vcf;
done < $1;
