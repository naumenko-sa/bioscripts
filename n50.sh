#!/bin/bash

#PBS -l walltime=100:00:00,nodes=1:ppn=12
#PBS -d .

if [ $# -ne "1" ]
then
    echo "Calculates N50 from contigs.fasta file"
    echo "Usage: n50.sh contigs.fasta"
    exit 1
fi

cat $1 | grep -v '>' | awk '{print length($0)}' | sort -rn > $1.len.sorted
echo "Avr len:" `cat $1.len.sorted | awk '{sum+=$0}END{print sum/NR}'`

len=`cat $1.len.sorted | awk '{sum+=$0}END{print sum}'`
echo "Total length: "$len
echo "N contigs: "`grep -c '>' $1`

len=$((len/2))
#echo $len

declare -i s=0;
cat $1.len.sorted | while read l;do s=$((s+l)); if [ $s -ge $len ]; then echo "N50: "$l; exit 0;fi; done; 

