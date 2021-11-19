#!/bin/bash

rm $1.tax
while read line;
do 
    id=`echo $line | awk '{print $2}' | awk -F '|' '{print $2}'`;
    echo -ne $line"\t" >> $1.tax;
    blast_get_taxonomy_gene.sh $id >> $1.tax;
done < $1;