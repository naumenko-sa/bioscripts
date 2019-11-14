#!/bin/bash

# https://slurm.schedmd.com/sbatch.html
# https://wiki.rc.hms.harvard.edu/display/O2

#SBATCH --partition=priority        # Partition (queue) priority
#SBATCH --time=2-00:00              # Runtime in D-HH:MM format, 10:00:00 for hours
#SBATCH --job-name=project          # Job name
#SBATCH -c 5			    # cores
#SBATCH --mem=50G           # Memory needed per CPU or --mem-per-cpu
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL             # Type of email notification (BEGIN, END, FAIL, ALL)

date

# render markdown
# conda activate r is not working on O2 nodes
# source activate r
# Rscript --vanilla -e 'rmarkdown::render("03.clustering_new.Rmd")'
# source deactivate

date