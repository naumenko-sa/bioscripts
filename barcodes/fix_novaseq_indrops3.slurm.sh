#!/bin/bash

# https://slurm.schedmd.com/sbatch.html

#SBATCH --partition=priority        # Partition (queue)
#SBATCH --time=2-00:00              # Runtime in D-HH:MM format
#SBATCH --job-name=bcbio       # Job name
#SBATCH -c 1			    # cores
#SBATCH --mem-per-cpu=15G           # Memory needed per CPU or --mem
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL             # Type of email notification (BEGIN, END, FAIL, ALL)

date

echo "file:" $1
echo "expand to": $2

gunzip -c $1 | awk -v full_length=$2 -f ~/code/bioscripts/barcodes/fix_novaseq_indrops3.awk | bgzip > $1.expanded.gz

date