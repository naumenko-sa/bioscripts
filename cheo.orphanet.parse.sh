#!/bin/bash

read_dom () {
    local IFS=\>
    read -d \< ENTITY CONTENT
}

while read_dom; 
do
    echo "$ENTITY => $CONTENT"
done < $1 > $1.txt

cat $1.txt | egrep "(Disorder id)|(Name lang)|ENSG" > $1.parsed

cat $1.parsed | awk '{if($0~"Disorder") {print $0;getline;print $0;}; if ($0~"ENSG") print $0;}' | grep -v "Disorder id" | awk -F '=> ' '{print $2}' > $1.final

cat $1.final | awk '{if($0 ~ "ENSG") {print $0"\t"dis}else{dis=$0}}' > orphanet.txt

#there are cases of many diseases for single gene
cat orphanet.txt | sort -k1,1 > orphanet.sorted.txt
cat orphanet.sorted.txt  | awk -F "\t" 'BEGIN{prev_gene="ZZ";buf=""}{if(prev_gene != $1){print prev_gene"\t"buf;buf=$2;prev_gene=$1}else{buf=buf","$2;}}END{print prev_gene","buf'} > orphanet.deduplicated.txt
