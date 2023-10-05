#!/bin/bash

#SBATCH --partition=priority        # Partition (queue) priority
#SBATCH --time=30-00:00             # Runtime in D-HH:MM format, 10:00:00 for hours
#SBATCH --job-name=project          # Job name
#SBATCH -c 1			    # cores
#SBATCH --mem=10G                  # total Memory or use --mem-per-cpu
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=NONE            # Type of email notification (BEGIN, END, FAIL, ALL)

date

tar czf $1.tar.gz $1
md5sum $1.tar.gz > $1.tar.gz.md5

date