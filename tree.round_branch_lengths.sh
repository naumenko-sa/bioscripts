#!/bin/bash

#rounds branch lengths in the tree file

cat $1 | sed -r s/":"/"\n"/g | sed 1d | awk -F ',' '{print $1}' | awk -F ')' '{print $1}' | awk '{printf $0"\t";printf("%.4f\n",$0);}' > $1.table
cp $1 $1.rounded

while read pattern rounded
do
    cat $1.rounded | awk -v p1=$pattern -v p2=$rounded '{gsub(p1,p2,$0);print $0;}' > $1.tmp;
    mv $1.tmp $1.rounded;
done < $1.table