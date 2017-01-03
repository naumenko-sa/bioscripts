#!/bin/bash

#list successfully finished bcbio runs only with family name
#for f in go.bcbio_array*;do cheo.check_if_done.sh $f >> ready.txt;done;
#while read log family;do mv $family ../0done/;done < ready.txt
#while read family;do qsub ~/bioscripts/cheo.postprocess.sh -v family=$family;done < families.txt 
#for f in `cat families.txt`;do mkdir 2report/$f;cp $f/*.txt 2report/$f;cp $f/*.table 2report/$f;done;

if grep -q "Timing: finished" $1;
then
    echo $1 `cat $1 | sed -n 3p | awk -F '/' '{print $8}'`
fi