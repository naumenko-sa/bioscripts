#!/bin/bash

# $1 - result of bam.coverage.sh with -d for a bed file
# calculates median coverage for a gene, and % of bases covered at 20x

echo "START:" `date`

cat $1 | awk '{print $4}' | sort | uniq  > genes.txt

echo "Parsing gene coverage:" `date`

for gene in `cat genes.txt`
do
    cat $1 | awk -v gen=$gene '{if($4 ~ gen ){print $0}} ' > $gene.genecov
done

echo "Calculating medians:" `date`

for f in *.genecov
do
    bam.coverage.median.sh $f
done

for gene in `cat genes.txt`
do
    cat $gene.genecov.median >> $1.medians
    rm $gene.genecov.median
done

echo "Calculating % of 20X coverage" `date`

for f in `cat genes.txt`
do 
    echo $f,`bam.coverage.base_wise.stat.py $f.genecov | grep "^20,"` | awk -F ',' '{print $4}'  >> $1.gene_wise
    rm $f.genecov
done

rm genes.txt

echo "END:" `date`
