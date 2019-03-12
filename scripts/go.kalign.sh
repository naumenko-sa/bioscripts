#!/bin/bash

kalign -gpo 500 -gpe 3 -tgpe 3 -bonus 0 -i $1 -o $1.fa -f fasta
align2fasta.pl $1.fa > $1.fasta
rm $1.fa