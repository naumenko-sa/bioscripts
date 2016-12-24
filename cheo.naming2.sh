#!/bin/bash

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

#process C4R names of type2:
#267626_GACTAGTA_R1_10.fastq.gz

#parameters project
echo "Project: "$project

target_path=/home/naumenko/work/project_c4r/data

original_path=/hpf/largeprojects/ccm/forgedata/c4r_unsolved_clinical_exome

cd ${original_path}/$project

for sample in `ls *R1*fastq.gz | awk -F '_' '{print $1}' | sort | uniq`
do
    echo "Sample: "$sample
    for pair in {R1,R2}
    do
	for f in `ls -v $sample*"_"$pair"_"*gz`
	do 
	    echo "File: "$f
	    gunzip -c $f >> $target_path/$project.${sample}_${pair}.fq
	done
	gzip $target_path/$project.${sample}_${pair}.fq
    done
done

cd -
