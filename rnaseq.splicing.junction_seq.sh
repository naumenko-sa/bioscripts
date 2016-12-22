#!/bin/bash

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=15g,mem=15g

module load R
Rscript ~/bioscripts/rnaseq.splicing.junction_seq.R