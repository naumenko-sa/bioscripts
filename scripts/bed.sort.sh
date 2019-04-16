#!/bin/bash

bedtools sort -faidx /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa.bedtoolsindex -i $1 > `echo $1 | sed s/.bed/.sorted.bed/`

