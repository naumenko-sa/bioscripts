#!/bin/bash

if [ $# -ne "5" ]
then
    echo "Runs trimmomatic"
    echo "Usage: trim_single.sh single.fq adapters.fasta quality[0-40] minlen threads"
    exit 1
fi

java -jar /home/tools/trimmomatic/trimmomatic-0.30.jar SE -threads $5 -phred33 \
  $1  \
  $1.trim.fq \
  ILLUMINACLIP:$2:2:40:15 LEADING:$3 TRAILING:$3 SLIDINGWINDOW:4:$3 MINLEN:$4
  

wc -l $1.trim.fq | awk '{print $1/4}' > $1.trim.fq.nreads
