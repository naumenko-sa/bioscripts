#!/bin/bash

# obtain a fresh copy of OMIM by email (request on the OMIM website)

# molecular basis for the disorder is known, mutations found in the gene
cat genemap2.txt | grep -v "^#" | grep '(3)' | grep ENSG | awk -F "\t" '{print $9"\t"$11"\t"$13}' | sort -k1,1 > omim.tmp

nlanger=`cat omim.tmp | grep -c ENSG00000185960`

if [[ $nlanger -gt 1 ]];then
    echo "Bug in OMIM: ENSG00000185960 still has 2 entries. Selecting one ... "
fi

# as of 2025, there is one gene with two idential entries but disease associations are in different order: ENSG00000185960
cat omim.tmp  | awk -F "\t" 'BEGIN{prev_gene="Ensembl_gene_id";buf="Omim_gene_description"}{if(prev_gene != $1){print prev_gene"\t"buf;buf=$2;prev_gene=$1}else{buf=buf","$2;}}END{print prev_gene"\t"buf'} > omim.txt

# 4 genes in omim (3) category don't have gene names: TRAC, IGHG2, IGHM, IGKC

#not generating omim_inheritance - needs manual curation
#cre.omim.inheritance.py genemap2.txt > omim_inheritance.txt

#cp omim.txt ~/cre
#cp omim_inheritance.txt omim.txt ~/cre

#rm omim.tmp