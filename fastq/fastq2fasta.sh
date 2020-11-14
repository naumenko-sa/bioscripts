#!/bin/bash

base=`echo $1 | sed s/.fq//`
echo "fq2fa:" $base
rm ${base}.fa
fastq2fasta -i $1 -o ${base}.fa
#rm ${base}.fasta.qual
