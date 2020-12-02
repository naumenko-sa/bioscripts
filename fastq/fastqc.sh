#!/bin/bash

#SBATCH --partition=priority        # Partition (queue) priority
#SBATCH --time=2-00:00              # Runtime in D-HH:MM format, 10:00:00 for hours
#SBATCH --job-name=qc               # Job name
#SBATCH -c 10			    # cores
#SBATCH --mem=20G                   # Memory needed per CPU or --mem-per-cpu
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=NONE             # Type of email notification (BEGIN, END, FAIL, ALL)


fastqc -noextract -nogroup -t 10 $1 $2 $3 $4
