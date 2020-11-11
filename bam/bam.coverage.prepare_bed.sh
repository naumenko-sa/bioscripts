#!/bin/bash

bedtools sort -faidx /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa.bedtoolsindex -i $1 | awk -F "\t" '{print $1"\t"$2-1"\t"$3"\t"$4}' > `echo $1 | sed s/.bed/.prepared.bed/`

