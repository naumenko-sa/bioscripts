#!/bin/bash
date
#remove reads w external adapters : 
adapters=( GATCGGAAGAGCACACGTCTGAACTCCAGTCAC GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT CTGTCTCTTATACACATCT AGATGTGTATAAGAGACAG );

wc -l $1 > before.wc

for a in "${adapters[@]}"
do

    mp_check_adapter.pl $1 $2 $a print 00

    rm $1 $2 

    mv $1.filt $1;
    mv $2.filt $2

done

wc -l $1 > after.wc
date