#!/bin/bash

# download: https://sites.google.com/site/revelgenomics/downloads
unzip revel-v1.3_all_chromosomes.zip
# tabulate and comment the header
cat revel_with_transcript_ids | tr "," "\t" | sed '1s/.*/#&/' > tabbed_revel.tsv
head -n1 tabbed_revel.tsv > revel_grch38.tsv
cat tabbed_revel.tsv | sed 1d | awk '$3 != "." ' | sort -k1,1 -k3,3n -T . >> revel_grch38.tsv
bgzip revel_grch38.tsv
tabix -f -s 1 -b 3 -e 3 revel_grch38.tsv.gz
rm revel_with_transcript_ids tabbed_revel.tsv
