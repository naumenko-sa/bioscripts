#!/bin/bash

# molecular basis for the disorder is known, mutations found in the gene
cat genemap2.txt | grep -v "^#" | grep '(3)' | grep ENSG | awk -F "\t" '{print $11"\t"$13}' | sort -k 1,1 > omim.tmp

#some genes may have two entries

cat omim.tmp  | awk -F "\t" 'BEGIN{prev_gene="Ensembl_gene_id\tOmim_gene_description";buf=""}{if(prev_gene != $1){print prev_gene"\t"buf;buf=$2;prev_gene=$1}else{buf=buf","$2;}}END{print prev_gene"\t"buf'} > omim.txt

#not generating omim_inheritance - needs manual curation
#cre.omim.inheritance.py genemap2.txt > omim_inheritance.txt

cp omim.txt ~/cre
#cp omim_inheritance.txt omim.txt ~/cre

rm omim.tmp