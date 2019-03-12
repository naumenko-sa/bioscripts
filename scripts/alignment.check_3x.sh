#!/bin/bash

#how many sequences in the alignment are not 3x in length

cat $1 | grep  -v '>' | awk '{print length($0)/3}' | grep '\.' | wc -l