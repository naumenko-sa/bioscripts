#!/bin/bash

# print 100 most frequent UMIs (first 8 bp) in $1

gunzip -c $1 | head -n 4000000 | awk '{if(NR%4==2) print substr($0,1,8)}' | sort | uniq -c |  awk '{print $2"\t"$1}' | sort -k2,2nr |  head -n100