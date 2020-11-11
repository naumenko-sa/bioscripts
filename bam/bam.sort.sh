#!/bin/bash -l

#SBATCH --job-name=samtools
#SBATCH --mem=8G
#SBATCH --export=ALL
#SBATCH -t 7-50:00
#SBATCH -p core -n 1

# sam2bam
# samtools view -bS $1 > $bamname

sortedbam=`echo $1 | sed s/bam/sorted.bam/`
samtools sort $bam $sortedbam
samtools index $sortedbam

