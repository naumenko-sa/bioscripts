#!/bin/bash -l

#SBATCH --job-name=crispresso
#SBATCH --mem=20G
#SBATCH --export=ALL
#SBATCH -t 7-50:00
#SBATCH -p core -n 10


date
# $1 = sample
# ml bcbio-nextgen/1.2.8

# otherwise r1 and r2 are not even
# bamtofastq filename=${1}-ready.bam F=${1}_1.fq F2=${1}_2.fq O=${1}_unmatched_1.fq O2=${1}_unmatched_2.fq S=${1}_single.fq
singularity run -e crispresso2_latest.sif CRISPRessoPooled -r1 ${1}_1.fq -r2 ${1}_2.fq -f baits.mm10.sorted.split.tsv -p 10 -n ${1}_small_baits

date
