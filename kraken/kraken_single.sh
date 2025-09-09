#!/bin/bash

# https://slurm.schedmd.com/sbatch.html
# https://wiki.rc.hms.harvard.edu/display/O2

#SBATCH --partition=priority        # Partition (queue) priority
#SBATCH --time=12:00:00              # Runtime in D-HH:MM format, 10:00:00 for hours
#SBATCH --job-name=kraken          # Job name
#SBATCH -c 10		           # cores
#SBATCH --mem=50G                  # Memory
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=NONE            # Type of email notification (BEGIN, END, FAIL, ALL)

# $1 - sample_name
# $2 - db path
# $3 - r1.fq.gz
# $4 - r2.fq.gz

/data/tools/kraken2bin/kraken2 \
--use-names \
--output - \
--threads 1 \
--report $1.single.kraken.txt \
--gzip-compressed \
--db $2 \
$3
#--classified_out $1.classified.txt \
#$2 $3
