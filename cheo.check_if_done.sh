#!/bin/bash

#list successfully finished bcbio runs only with family name

if grep -q "Timing: finished" $1;
then
    echo $1 `cat $1 | sed -n 3p | awk -F '/' '{print $8}'`
fi