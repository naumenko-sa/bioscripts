#!/bin/bash

# $1 = input.fasta
# output = input.internal.fasta

prefix=~/project_reversals/scripts/paml
filename=`echo $1 | sed s/fasta/phy/`
perl ~/bin/fasta2phylip.pl $1 $filename
cat $prefix/codeml.template | sed s/sequence.phy/$filename/ > codeml.ctl
codeml
pamlresult=`echo $1 | sed s/fasta/paml/`
mv rst $pamlresult
#rm 2base.t in.basemlg rates rst1 rub baseml.ctl lnf mlb $filename;
perl $prefix/paml_parse.pl $pamlresult > `echo $1 | sed s/fasta/internal.fasta/`
