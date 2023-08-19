#!/bin/bash

#SBATCH -t 7-00:00:00
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --job-name=cram2bam
#SBATCH --mail-type=NONE

date

cramtools \
bam \
-R /n/holylfs05/LABS/hsph_bioinfo/Everyone/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa \
-I $1 \
--output-bam-file $2 \
-b

date
