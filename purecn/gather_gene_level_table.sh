#!/bin/bash

# generate a table of gene-level CNA from bcbio purecn output

#cd final/project

cat */purecn/*_genes.csv | head -n1
cat */purecn/*_genes.csv | grep -v Sampleid | awk -F ',' 'BEGIN{OFS=","}{if ($6!=2) printf "%s,%s,%s,%s,%s,%s,%s,%.2f,%s,%s,%.2f,%.2f,%.2f,%s,%s,%s,%s,%s,%s,%s\n", $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20}' 
