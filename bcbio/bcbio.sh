#!/bin/bash

# https://slurm.schedmd.com/sbatch.html

#SBATCH --partition=priority        # Partition (queue)
#SBATCH --time=10-00:00:00             # Runtime in D-HH:MM format
#SBATCH --job-name=project         # Job name
#SBATCH -c 20
#SBATCH --mem-per-cpu=5G            # Memory needed per CPU
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL             # Type of email notification (BEGIN, END, FAIL, ALL)

bcbio_nextgen.py ../config/project.yaml -n 20
