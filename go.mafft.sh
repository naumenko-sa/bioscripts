#!/bin/bash

mafft --quiet --thread 20 --preservecase --auto $1 > $1.aln;
align2fasta.pl $1.aln > `echo $1 | sed s/fasta/aln.fasta/`;
rm $1.aln;