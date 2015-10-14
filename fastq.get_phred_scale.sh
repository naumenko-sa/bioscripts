#!/bin/bash

head -n 4000 $1 | awk '{if (NR%4==0) print $0;}' | perl -e 'while(<>){for(my $i=0;$i<length($_)-1;$i++){print ord(substr($_,$i,1))."\n"};}' | sort -rn > $1.values
echo "Min phred:".`tail -n1 $1.values`;
echo "Max phred:".`head -n1 $1.values`;