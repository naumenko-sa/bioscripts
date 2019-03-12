#!/bin/bash

#2016-09-30: filter GEMINI output for excel file

bname=`echo $1 | sed s/.txt//`

cat $1 | egrep -v "^chrom|RYR1|CACNA1S|ASPH" > ${bname}.othergenes;

#if the variant is rare and exonic,or splicing
head -n1 $bname.othergenes > $bname.filtered.tmp1;
#field 13 = is_exonic
#field 16 = is_splicing
cat ${bname}.othergenes | awk '{if ($36<0.01 && ($13==1 || $16==1)) print $0;}' >> $bname.filtered.tmp1;

#cat $1.3genes | egrep "possibly_damaging|probably_damaging" >>$1.filtered.tmp;
#			 head -n1 $bname.filtered.tmp | cut -f 1,3,5,6,10,11,12,17,19,20,24-27,29,33,36,37 > $bname.filtered;
#cat $bname.filtered.tmp | grep -v chrom | sort | uniq | cut -f 1,3,5,6,10,11,12,17,19,20,24-27,29,33,36,37 >> $bname.filtered;
#rm $bname.filtered.tmp;