#!/bin/bash

#PBS -l walltime=1000:00:00,nodes=1:ppn=1
#PBS -d .

build_lmer_table  -sequence contig_150_449.nononplant.fasta -freq freq -v
repeatscout -sequence contig_150_449.nononplant.fasta -output repeats.fasta -freq freq -vvvv
cat repeats.fasta | filter-stage-1.prl > repeats.filtered.fasta
repeatmasker -lib repeats.filtered.fasta  -gff contig_150_449.nononplant.fasta