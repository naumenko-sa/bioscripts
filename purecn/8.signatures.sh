#!/bin/bash -l

#SBATCH --partition=short          # Partition (queue) priority
#SBATCH --time=10:00:00            # Runtime in D-HH:MM format, 10:00:00 for hours
#SBATCH --job-name=purecn          # Job name
#SBATCH -c 1                       # cores
#SBATCH --mem=10G                  # Memory needed per CPU or --mem-per-cpu
#SBATCH --output=project_%j.out    # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err     # File to which STDERR will be written, including job ID
#SBATCH --mail-type=NONE           # Type of email notification (BEGIN, END, FAIL, ALL)


# #SBATCH --job-name=purecn
# #SBATCH --mem=30G
# #SBATCH --export=ALL
# #SBATCH -t 7-50:00
# #SBATCH -p core -n 10

# once OOM killed with 20G and bootstrap500
date

. .profile
which Rscript

# $1 = sample-callable.bed
# $2 = sample_coverage_loess.RDS

SAMPLEID=`echo $2 | awk -F '_' '{print $1}'`

echo $SAMPLEID

Rscript $PURECN/Dx.R \
--out $SAMPLEID \
--rds $2 \
--callable $1 \
--exclude $PURECN/hg38_simpleRepeats.bed \
--signatures

date
