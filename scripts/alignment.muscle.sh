#!/bin/bash

# align a fasta file with muscle
# http://www.drive5.com/muscle/manual/

muscle -in $1 -out _tmp.fa
alignment.fa2fasta.pl _tmp.fa > $1.aln
rm _tmp.fa

python ~/bioscripts/scripts/muscle.stabilize.py $1 $1.aln > $1.stab

rm $1.aln
mv $1.stab `echo $1 | sed s/fasta/align.fasta/`

