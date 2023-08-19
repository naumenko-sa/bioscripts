#!/bin/bash

# max 3 days
#SBATCH -t 3-00:00:00
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --job-name=cram2fq
#SBATCH --mail-type=NONE

date

cramtools \
fastq -R /n/holylfs05/LABS/hsph_bioinfo/Everyone/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa \
--gzip \
-I $1 \
--fastq-base-name $2

date
