#!/bin/bash

#calculates R/A(1)/R/A(3)statistics

kstat()
{
    a=`cat $1 | awk -v p=" $2$" '{if ($0 ~ p){getline;getline;print $5;getline;getline;print $5;}}'| head -n1`;
    b=`cat $1 | awk -v p=" $2$" '{if ($0 ~ p){getline;getline;print $5;getline;getline;print $5;}}'| tail -n1`; 

    echo $a' '$b | awk '{if ($2 != 0){printf "%.3f\n",$1/$2} else {print $2}}';
}

#$1=10mln.1000.0.001.1;
rm $1.stat;
for f in `seq 10000 10000 $2`;
do
    kstat $1 $f >> $1.stat;
done;

