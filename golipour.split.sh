#!/bin/bash

#cat $1 | awk '{if (NR%2==1){ rname=$1;}else {print rname" /1";print substr($0,1,100);print rname" /2";print substr($0,103,100); }}' > temp.fq;

cat temp.fq | awk '{if ($0 ~ "\/1$"){print $0;getline;print $0;}}' > left.fq;

cat temp.fq | awk '{if ($0 ~ "\/2$"){print $0;getline;print $0;}}' > right.fq;

rm temp.fq;