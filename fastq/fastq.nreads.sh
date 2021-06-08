#!/bin/bash

#SBATCH --partition=short           # Partition (queue) priority
#SBATCH --time=12:00              # Runtime in D-HH:MM format, 10:00:00 for hours
#SBATCH --job-name=qc               # Job name
#SBATCH -c 1			    # cores
#SBATCH --mem=10G                   # Memory needed per CPU or --mem-per-cpu
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=NONE             # Type of email notification (BEGIN, END, FAIL, ALL)


zcat $1 | wc -l | awk '{print $0/4}' > $1.nreads
