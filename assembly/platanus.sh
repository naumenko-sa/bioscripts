#!/bin/bash

# http://platanus.bio.titech.ac.jp/

#PBS -d .
#PBS -l walltime=1000:00:00,nodes=node01:ppn=48,mem=520GB

date

#platanus assemble -t 24  -f *.fastq 2> ass_log.txt
platanus scaffold -c $1  -op1 $2  &> scaff_log.txt
platanus gap_close -c out_scaffold.fa -op1 $2 &> gap_close.txt

#platanus bubble_map -c out_gapClosed.fa -b out_contigBubble.fa out_scaffoldBubble.fa 2>bub_log.txt

date
