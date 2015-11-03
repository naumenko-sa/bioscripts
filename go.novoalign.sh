#!/bin/bash

base=/mnt/lustre/tools/novocraft
#${base}/novoalign -d hla.clean -f $1 $2 -t 30 -o SAM -r All 5 -l 80 -e 100 -i PE 600 200  > aln.sam
${base}/novoalign -d contig -f $1 $2 -t 150 -e 100 -r All 10 -l 50 -o Pairwise -v 300 -q 2 -i PE 50-600 > aln.sam
#sam2sortedbam.sh aln.sam
