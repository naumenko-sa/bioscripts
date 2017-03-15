#!/bin/bash

# reverse translates amino acid alignment to DNA
# uses faSomeRecords from Kent tools and revtrans
# cds.fasta is a database of all cds sequences

# $1 - name.whatever.fasta, i.e. name.align.fasta, protein alignment
# cds.fasta in the current dir = database of all cds sequences
# result = name.fasta

bname=`echo $1 | awk -F "." '{print $1}'`

cat $1 | grep '>' | sed s/">"// > $bname.names

#makes unaligned file with CDS sequences
faSomeRecords cds.fasta $bname.names $bname.cds

/home/tools/RevTrans-1.4/revtrans.py -match pos $bname.cds $1 $bname.fa;

alignment.fa2fasta.pl $bname.fa > $bname.fasta

rm $bname.fa $bname.names $bname.cds
