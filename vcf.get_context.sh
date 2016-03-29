#!/bin/bash

#for bhutani c-0
#for dbsnp c-1

while read line;
do
    read chr pos ref alt <<< "$line";
    triplet=`tail -n1  chr${chr}.fasta | awk -v c=$pos '{print substr($0,c-1,3)}'`;
    echo $chr $pos $ref $alt $triplet;
done < $1;