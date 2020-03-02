#!/bin/bash

# https://slurm.schedmd.com/sbatch.html

#SBATCH --partition=short        # Partition (queue)
#SBATCH --time=02:00:00              # Runtime in D-HH:MM format
#SBATCH --job-name=fastqc            # Job name
#SBATCH -c 10
#SBATCH --mem-per-cpu=1G            # Memory needed per CPU
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL             # Type of email notification (BEGIN, END, FAIL, ALL)

which java
which fastqc
fastqc --noextract --nogroup -t 10 $1
