#!/bin/bash

if [ $# -lt 4 ];
then
    echo "Usage: `basename $0` r1.fq r2.fq minqual minlength [adapter]";
    exit;
fi;

adapt='';

if [ $# -eq 5 ];#adapter seq
then
    adapt="--pcr1 $5 --pcr2 $5";
fi;

base=`echo $1 | sed s/.r1.fq//`;

date
adapterremoval --file1 $1 --file2 $2 --basename $base  --trimns --trimqualities --minquality $3 --collapse --minlength $4 $adapt --stats
rm $base.discarded $base.collapsed
wc -l $base.collapsed.truncated > $base.collapsed.truncated.wc
wc -l $base.pair1.truncated > $base.pair1.truncated.wc
wc -l $base.singleton.truncated > $base.singleton.truncated.wc

echo "Paired: "`cat $base.pair1.truncated.wc`;
echo "Collapsed: "`cat $base.collapsed.truncated.wc`;
echo "Singleton: "`cat $base.singleton.truncated.wc`;
date