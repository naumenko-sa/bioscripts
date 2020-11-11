#!/bin/bash -l

#SBATCH --job-name=samtools
#SBATCH --mem=8G
#SBATCH --export=ALL
#SBATCH -t 7-50:00
#SBATCH -p core -n 1

# sam2bam
# samtools view -bS $1 > $bamname

# $1 = sample.unsorted.bam

bname=`basename $1 .bam`

samtools sort $1 $bname.sorted.bam
samtools index $bname.sorted.bam

