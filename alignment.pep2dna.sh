#!/bin/bash

# reverse translates amino acid alignment in to DNA
# uses faSomeRecords from Kent tools and revtrans
# cds.fasta is a database of all cds sequences

# $1 - name.aln.fasta
# result = name.cds.aln

cat $1 | grep '>' | sed s/">"// > $1.names

cds_file=`echo $1 | sed s/aln/cds/`

#makes unaligned file with CDS sequences
faSomeRecords cds.fasta $1.names $cds_file

/home/tools/RevTrans-1.4/revtrans.py -match pos $cds_file $1 $1.fa;

align2fasta.pl $1.fa >`echo $1 | sed s/aln/cds.aln/`;

rm $1.fa $1.names $cds_file