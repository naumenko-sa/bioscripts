#!/bin/bash

# https://slurm.schedmd.com/sbatch.html

#SBATCH --partition=priority        # Partition (queue)
#SBATCH --time=5-00:00:00           # Runtime in D-HH:MM format
#SBATCH --job-name=bcbio            # Job name
#SBATCH -c 1
#SBATCH --mem-per-cpu=10G           # Memory needed per CPU
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=NONE            # Type of email notification (BEGIN, END, FAIL, ALL)

bcbio_nextgen.py ../config/bcbio.yaml -n 72 -t ipython -s slurm -q medium -r t=0-72:00 -r conmem=20 --timeout 3000 --tag "rna"
