#!/bin/bash

#makeblastdb -in UniVec_Core -dbtype 'nucl'
go.blastn.sh $1 contaminants.fasta 100 1e-1
cat ${1}_vs_contaminants.fasta.blastn.1e-1 | awk '{print $1}' | sort | uniq > $1.contaminated
faSomeRecords -exclude $1 $1.contaminated $1.clean

go.blastn.sh $1.clean UniVec 100 1e-1
cat $1.clean_vs_UniVec.blastn.1e-1 | awk '{print $1}' | sort | uniq > $1.clean.contaminated
faSomeRecords -exclude $1.clean $1.clean.contaminated $1.clean.clean

rm $1.clean
mv $1.clean.clean $1.clean
rm $1.contaminated $1.clean.contaminated *.blastn.1e-1
