#!/bin/bash -l

#SBATCH --job-name=preseq
#SBATCH --mem=8G
#SBATCH --export=ALL
#SBATCH -t 7-50:00
#SBATCH -p core -n 1


# http://smithlabresearch.org/software/preseq/

# $1 = file.bam

bname=`basename $1 .bam`

preseq c_curve -bam -o $bname.out $1
