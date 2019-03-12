#!/bin/bash

median_line=`cat $1 | wc -l`
median_line=$(($median_line/2))
cat $1 | awk '{print $6}' | sort -n | sed -n ${median_line}p > $1.median
