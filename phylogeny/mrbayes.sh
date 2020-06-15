#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=100:00:00,mem=10gb
#PBS -d.
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/mnt/lustre/tools/beagle-lib/lib
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64/openmpi/lib
#mpirun -np 8 mb test.nex

/home/tools/mrbayes/3.2.6/bin/mb $nex