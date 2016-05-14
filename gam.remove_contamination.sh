#!/bin/bash

#makeblastdb -in UniVec_Core -dbtype 'nucl'
go.blastn.sh $1 $2 100 1e-1
cat ${1}_vs_${2}.blastn.1e-1 | awk '{print $1}' | sort | uniq > $1.contaminated
faSomeRecords -exclude $1 $1.contaminated $1.clean