#!/bin/bash

#just compare 2 vcf files - no vcf output

module load vcftools

vcftools --vcf $1 --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 \
         --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 --chr X --chr Y --chr M --diff $2 --diff-site --out $1.$2
