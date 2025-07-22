#!/bin/bash

# max 3 days
#SBATCH -t 3-00:00:00
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --job-name=cram2fq
#SBATCH --mail-type=NONE

date

# $1 sample.cram
# $2 sample_name
# $3 reference.fa

cramtools \
fastq -R $3 \
--gzip \
-I $1 \
--fastq-base-name $2

date
