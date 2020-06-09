#!/bin/bash
#PBS -l walltime=10:00:00,nodes=1:ppn=5
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

# modeltest.sh file.fasta num_threads
java -jar /n/data1/cores/bcbio/naumenko/tools/jmodeltest-2.1.10/jModelTest.jar -d $1 -s 3 -g 4 -f -AIC -BIC -i -a -tr 5 -o $1.modeltest
#-s 11 -i
