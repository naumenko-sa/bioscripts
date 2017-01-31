#!/bin/bash

#prepares a family to be merged to SEEN_IN_C4R database

if [ -z $family ]
then
    family=$1
fi

cd $family

nsamples=`cat samples.txt | wc -l`

nfield=5
for sample in `cat samples.txt`;
do
    cat $family.txt | cut -d ";" -f 1,3,4,$nfield | grep -v ";-" | grep -v "Insufficient" | awk -F ";" -v sam=$sample '{print $1"-"$2"-"$3"\t"sam}' | sed 1d | sort -k1,1 > $family.$sample.c4r
    nfield=$(($nfield+1))
done

cd ..
