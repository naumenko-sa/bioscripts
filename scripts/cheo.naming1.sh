#!/bin/bash

#process filenames of C4R of type1: [sample].pair1.fastq.gz

#parameters $1=project1

target_path=/home/naumenko/work/project_c4r/data

original_path=/hpf/largeprojects/ccm/forgedata/c4r_unsolved_clinical_exome

cd ${original_path}/$1

for f in *.pair1.fastq.gz
do
    sample=`echo $f | awk -F "." '{print $1}'`
    ln -s $original_path/$1/$f $target_path/$1.${sample}_1.fq.gz
    fastq2=`echo $f | sed s/pair1/pair2/`
    ln -s $original_path/$1/$fastq2 $target_path/$1.${sample}_2.fq.gz
done

cd -
