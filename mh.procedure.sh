#!/bin/bash

#filter variants for MH
#3 genes
gunzip -c $1/final/$1/$1-gatk-haplotype.vcf.gz | egrep "RYR1|CACNA1S|ASPH" > $1.3genes
vcfvep2txt.sh $1.3genes
#canonical isoform
cat $1.3genes.txt | egrep "YES|CHR" > $1.3genes.canonical
#3 genes
cat $1/final/2016-08-30_mh/$1-gatk-haplotype.txt | egrep "RYR1|CACNA1S|ASPH" > $1.gem
#exac frequency
cheo.gemini_vep_merge.pl $1.3genes.canonical $1.gem > $1.merged

cat $1.merged | egrep "CHR|RYR1|CACNA1S|ASPH" |  egrep "CHR|missense" | awk -F "\t" '{if ($23<0.01 || $1=="CHR")  print $0;}' > $1.tmp;

cat $1.tmp | awk -F "\t" '{if ($1=="CHR") print $0;else{split($9,a,"/");split($16,pos,"/");for(i=1;i<9;i++)printf $i"\t";printf a[1]pos[1]a[2];for (i=10;i<=NF;i++)printf "\t"$i;printf "\n";}}' > $1.txt;

rm $1.3genes $1.gem $1.3genes.txt $1.3genes.canonical  $1.merged $1.tmp

