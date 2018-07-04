#!/bin/bash

#$1 - result of bam.coverage.sh with -d for a bed file

cat $1 | awk '{print $4}' > genes.txt

for gene in `cat genes.txt`
do 
    $1 | awk -v gen=$gene '{if($4 ~ gen ){print $0}} ' > $gene.genecov;
done

for f in *.genecov
do 
    bam.coverage.median.sh $f
done

for f in `cat genes.txt | sort`
do 
    echo $f,`bam.coverage.base_wise.stat.py $f.genecov | grep "^20,"` >> $1.gene_wise
done

