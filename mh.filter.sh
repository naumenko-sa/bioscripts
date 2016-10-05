#!/bin/bash

#2016-09-30: filter GEMINI output for excel file

cat $1 | egrep "^chrom|RYR1|CACNA1S|ASPH" > $1.3genes;
#if the variant is rare and exonic
head -n1 $1.3genes > $1.filtered;
cat $1.3genes | awk '{if ($36<0.01 && $13==1) print $0;}' > $1.filtered.tmp;
#cat $1.3genes | egrep "possibly_damaging|probably_damaging" >>$1.filtered.tmp;
cat $1.filtered.tmp | sort | uniq | cut -f 1,3,5,6,10,11,12,17,19,20,24-27,29,33,36,37> $1.filtered;
rm $1.filtered.tmp;