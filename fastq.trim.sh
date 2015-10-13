#!/bin/bash

if [ $# -ne "6" ]
then
    echo "Runs trimmomatic"
    echo "Usage: trim.sh left.fq right.fq adapters.fasta quality[0-40] minlen threads"
    exit 1
fi

java -jar /home/tools/trimmomatic/trimmomatic-0.30.jar PE -threads $6 -phred33 \
  $1 $2 \
  `echo $1 | sed s/fq/trim.fq/` forward_unpaired.fq `echo $2 | sed s/fq/trim.fq/` reverse_unpaired.fq \
  ILLUMINACLIP:$3:2:40:15 LEADING:$4 TRAILING:$4 SLIDINGWINDOW:4:$4 MINLEN:$5

nreads.sh `echo $1 | sed s/fq/trim.fq/`

cat forward_unpaired.fq reverse_unpaired.fq > unpaired.fq
nreads.sh unpaired.fq

rm forward_unpaired.fq reverse_unpaired.fq