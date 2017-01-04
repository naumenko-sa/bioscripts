#!/bin/bash

gunzip -c QC.spliceJunctionCounts.novelSplices.txt.gz | sed 1d | cut -f 1,3,4,5 > $1.novel_junctions.bed
bedtools intersect -wa -wb -a $1.novel_junctions.bed -b ~/Desktop/reference_tables/muscular_genes.bed | awk '{print $1"\t"$2"\t"$3"\t"$8"\t"$4}' | awk '{if ($5>5) print $0}' > $1.results.bed

