#!/bin/bash

/home/tools/trinityrnaseq_r20131110/trinity-plugins/transdecoder/TransDecoder -t $1
rm -rf transdecoder.tmp.*
rm $1.transdecoder.bed
rm $1.transdecoder.gff3

align2fasta.pl $1.transdecoder.pep | awk '{name=$0;name1=$1;name2=$6;getline;print name1" "name2"_"length($0)*3;print $0;}' | sed s/cds.// | sed s/\|.*type:/_/ > $1.pep
align2fasta.pl $1.transdecoder.cds | awk '{name=$0;name1=$1;name2=$6;getline;print name1" "name2"_"length($0);print $0;}' | sed s/cds.// | sed s/\|.*type:/_/ > $1.cds

rm $1.transdecoder.pep $1.transdecoder.cds 