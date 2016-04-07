#!/bin/bash

# check -10 and +10 conservative nucleotide positions in codon alignement $1 for amino acid position (1-based) $2
len=`cat $1 | grep -v '>' | head -n1 | awk '{print length/3 - 4}'`


result=0
if (( $2 < 5 )) || (( $2 > len));then
    result=0;
else
    #unique_strings=
    cat $1  | grep -v '>' | awk -v pos=$2 '{print substr($0,3*(pos-1)+1-10,23)}' | awk '{print substr($0,1,10)substr($0,14,10)}' 
    #| sort | uniq | wc -l`;
    if [ $unique_strings -eq 1 ]; then
	result=1;
    fi
fi

echo $result
#| sort | uniq | wc -l