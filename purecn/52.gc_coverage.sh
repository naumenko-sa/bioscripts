#!/bin/bash -l

#SBATCH --job-name=purecn
#SBATCH --mem=20G
#SBATCH --export=ALL
#SBATCH -t 7-50:00
#SBATCH -p core -n 1

date

. .profile
### 5.2 gc_coverage
# $1 = sample.bam
# $2 = intervals.txt = output of process_intervals

which Rscript

Rscript \
$PURECN/Coverage.R \
--outdir . \
--bam $1 \
--intervals $2

date
