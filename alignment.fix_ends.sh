#!/bin/bash 

#looks for aa substitution in first and last 16 codons in an alignment
#if there are >=2 substitutions those regions are filtered out
#requires transeq from EMBOSS and alignment.nucleotide_diversity.pl from bioscripts

start=1

len=`cat $1 | head -n2 | tail -n1 | awk '{print length}'`

cat $1 | awk '{print $0;getline; print substr($0,1,48);}' > $1.tmp

transeq  -sequence $1.tmp -outseq $1.tmp.aa

div=`alignment.nucleotide_diversity.pl $1.tmp.aa | awk '{sum+=$2}END{print sum}'`;

#in two columns there are more than 2 aa's
if (( $div > 17));
then
    start=49;
    len=$(($len-48));
    #echo $len
fi;

rm $1.tmp $1.tmp.aa

cat $1 | awk '{print $0;getline;print substr($0,length($0)-47,48);}' > $1.tmp
transeq -sequence $1.tmp -outseq $1.tmp.aa

div=`alignment.nucleotide_diversity.pl $1.tmp.aa | awk '{sum+=$2}END{print sum}'`

if (( $div > 17));
then
    len=$(($len-48));
    #echo $len
fi;

rm $1.tmp $1.tmp.aa

cat $1 | awk -v st=$start -v len=$len '{print $0;getline;print substr($0,st,len);}'

