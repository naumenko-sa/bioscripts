#!/bin/bash

# convert EMG exported table to VEP input
cat EMG\ output\ Sept\ 2024.csv | cut -d, -f 1-4,6 | grep -v "<p." | grep -v "<div" | grep "NM_" | sed s/","/"\t"/g > variants.tsv
# produce tab output (not good for indels)
#cat variants.tsv | awk -F '\t' '{print $1"\t"$2"\t"$2"\t"$3"/"$4"\t+\tNSOVAR"NR}' > variants.vep_input.txt

# produce vcf output - ok for indels
cat variants.tsv | awk -F '\t' '{print $1"\t"$2"\tNSOVAR"NR"\t"$3"\t"$4"\t.\t.\t."}' > variants.vep_input_vcf.txt
