#!/bin/bash

a=`cat $1 | awk -v p=" $2$" '{if ($0 ~ p){getline;getline;print $6;getline;getline;print $6;}}'| head -n1`;
b=`cat $1 | awk -v p=" $2$" '{if ($0 ~ p){getline;getline;print $6;getline;getline;print $6;}}'| tail -n1`; 

echo $a' '$b | awk '{printf "%.3f\n",$1/$2}';