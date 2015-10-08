#!/bin/bash

cat $1 | grep '>' | sed s/">"// > $1.names;
faSomeRecords cds.fasta $1.names `echo $1 | sed s/aln/cds/`;
/home/tools/RevTrans-1.4/revtrans.py -match pos `echo $1 | sed s/aln/cds/` $1 $1.fa;
align2fasta.pl $1.fa >`echo $1 | sed s/aln/cds.aln/`;
rm $1.fa $1.names `echo $1 | sed s/aln/cds/`;