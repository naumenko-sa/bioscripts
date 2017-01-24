#!/bin/bash

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=20g,mem=20g

module load java

cramtools fastq -Xmx10g -F $sample --skip-md5-check \
		-z \
		-I $cram \
		-R /hpf/largeprojects/ccm/forgedata/c4r_unsolved_clinical_exome/1092R/hg19.fa.gz
