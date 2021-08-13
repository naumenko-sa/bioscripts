#!/bin/bash

cat E*_T/purecn/E*_T_genes.csv | head -n1  | awk -F ',' '{print $1"\t"$2"\t"$6"\t"$(NF-4)"\t"$(NF-2)}' > ERBB2.tsv
cat E*_T/purecn/E*_T_genes.csv | egrep "ERBB2," | awk -F ',' '{print $1"\t"$2"\t"$6"\t"$(NF-4)"\t"$(NF-2)}' >>ERBB2.tsv
cat E*_T/purecn/E*_T_genes.csv | head -n1  | awk -F ',' '{print $1"\t"$2"\t"$6"\t"$(NF-4)"\t"$(NF-2)}' > PTEN.tsv
cat E*_T/purecn/E*_T_genes.csv | egrep "PTEN," | awk -F ',' '{print $1"\t"$2"\t"$6"\t"$(NF-4)"\t"$(NF-2)}' >>PTEN.tsv
