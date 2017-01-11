#!/bin/bash

#parses qorts output looking for novel junctions in muscular genes and novel exons like found in DMD2016 article
#$1 = sample name: muscle2

gunzip -c QC.spliceJunctionCounts.novelSplices.txt.gz | sed 1d | cut -f 1,3,4,5 > $1.novel_junctions.bed
bedtools intersect -wa -wb -a $1.novel_junctions.bed -b ~/Desktop/reference_tables/muscular_genes.bed | awk '{print $1"\t"$2"\t"$3"\t"$8"\t"$4}' | awk '{if ($5>5) print $0}' > $1.results.bed
#get novel exons $2 = 100 or 200
cat $1.results.bed | awk '{prev=$3;getline;print $0"\t"$2-prev}' | awk -v threshold=$2 '{if ($6>0 && $6<threshold) print $0;}' > $1.novel_exons

gunzip -c QC.spliceJunctionCounts.knownSplices.txt.gz | sed 1d | cut -f 2,4,5,6 > $1.known_junctions.bed
bedtools intersect -wa -wb -a $1.known_junctions.bed -b ~/Desktop/reference_tables/muscular_genes.bed | awk '{print $1"\t"$2"\t"$3"\t"$8"\t"$4}' | awk '{if ($5<2) print $0}' > $1.known.low_coverage