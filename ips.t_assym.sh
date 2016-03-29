#!/bin/bash

while read line;
do
    read chr t pos ref alt t t t <<<$line;
    result=`cat hg19.refgene.strand |  awk -v chr="$chr" -v pos="$pos" '{if($1 == chr && $2<=pos && pos<=$3) print $0;}'`;
    echo $chr $pos $ref $alt $result;
done < $1;