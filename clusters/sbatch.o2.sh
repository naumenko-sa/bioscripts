#!/bin/bash

# https://slurm.schedmd.com/sbatch.html
# https://wiki.rc.hms.harvard.edu/display/O2

#SBATCH --partition=priority        # Partition (queue) priority
#SBATCH --time=30-00:00              # Runtime in D-HH:MM format, 10:00:00 for hours
#SBATCH --job-name=project          # Job name
#SBATCH -c 20			    # cores
#SBATCH --mem=100G                  # total Memory or use --mem-per-cpu
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL             # Type of email notification (BEGIN, END, FAIL, ALL)

date

# render markdown
# conda activate r is not working on O2 nodes
# source activate r
# Rscript --vanilla -e 'rmarkdown::render("03.clustering_new.Rmd")'
# source deactivate

# memory usage of a running job
# sstat job_id --format AveRSS,MaxRSS

# efficienty of a finished job
# seff job_id

bcbio_nextgen.py ../config/bcbio.yaml -n 20

date