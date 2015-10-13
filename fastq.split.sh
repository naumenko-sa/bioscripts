#!/bin/bash

if [ $# -ne "3" ]
then
    echo "Splits fastq file into left and right"
    echo "Usage: splitfastq.sh whole.fq left.fq right.fq"
    exit 1
fi

cat $1 | awk '{if ((NR%8==1) || (NR%8==2) || (NR%8==3) || (NR%8==4)) print $0}' > $2;
cat $1 | awk '{if ((NR%8==5) || (NR%8==6) || (NR%8==7) || (NR%8==0)) print $0}' > $3;