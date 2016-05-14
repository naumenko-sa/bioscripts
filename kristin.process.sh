#!/bin/bash

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d /home/naumenko/work
#PBS -l vmem=20g

module load p7zip/15.14.1

FDIR=/hpf/largeprojects/ccmbio/boycott_exome/raw

echo "Start " `date`

mkdir /scratch/${file}_dir
cp $FDIR/$file /scratch/${file}_dir

echo "Copy done " `date`

cd /scratch/${file}_dir
7za e -bd $file


echo "Unpacking done" `date`

for f in *.gz;do gunzip $f;done;

cat *R1*> $file.r1.fq
cat *R2*> $file.r2.fq

gzip $file.r1.fq
gzip $file.r2.fq

echo "gzip done " `date`
mv $file.r1.fq.gz ~/work/kristin_agilent/input
mv $file.r2.fq.gz ~/work/kristin_agilent/input

rm /scratch/${file}_dir/*

echo "End " `date`