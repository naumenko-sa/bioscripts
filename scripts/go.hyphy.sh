#!/bin/bash

bname=`echo $1 | awk -F "." '{print $1}'`;
alignment.remove_stops_n_gaps.pl $1 | awk -F '|' '{print $1}' > align.fasta;
HYPHYMP <go.hyphy.sh.input | tail -n7 > $bname.hyphy
