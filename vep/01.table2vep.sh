#!/bin/bash

# convert EMG exported table to VEP input
# the first line also removed the header 
cat EMG\ output\ Sept\ 2024.csv | cut -d, -f 1-4,6 | grep -v "<p." | grep -v "<div" | grep "NM_" | sed s/","/"\t"/g > variants.tsv
cat variants.tsv | awk -F '\t' '{print $1"\t"$2"\t"$2"\t"$3"/"$4"\t+\tNSOVAR"NR}' > variants.vep_input.txt
