#!/bin/bash

# https://slurm.schedmd.com/sbatch.html
# https://wiki.rc.hms.harvard.edu/display/O2

#SBATCH --partition=priority         # Partition (queue) priority
#SBATCH --time=5-00:00              # Runtime in D-HH:MM format, 10:00:00 for hours
#SBATCH --job-name=tar          # Job name
#SBATCH -c 1			    # cores
#SBATCH --mem=20G                  # total Memory or use --mem-per-cpu
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=NONE             # Type of email notification (BEGIN, END, FAIL, ALL)

# $1 = dir_to_archive

date
#tar czf $1.tar.gz $1
#md5sum $1.tar.gz > $1.tar.gz.md5
#date

tar cf $1.tar $1
md5sum $1.tar > $1.tar.md5

date
