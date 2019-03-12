#!/bin/bash

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d /home/naumenko/work
#PBS -l vmem=20g

#parameters: from target

echo "Start " `date`

mkdir /scratch/kristin_dir
cp ${from}/* /scratch/kristin_dir

echo "Copy done " `date`

cd /scratch/kristin_dir

for f in *.gz;do gunzip $f;done;

cat *R1*> ${target}_1.fq
cat *R2*> ${target}_2.fq

gzip ${target}_1.fq
gzip ${target}_2.fq

echo "gzip done " `date`

mv ${target}_1.fq.gz  /hpf/largeprojects/ccmbio/cheo.variants/data/Illumina
mv ${target}_2.fq.gz  /hpf/largeprojects/ccmbio/cheo.variants/data/Illumina

rm /scratch/kristin_dir/*

echo "End " `date`