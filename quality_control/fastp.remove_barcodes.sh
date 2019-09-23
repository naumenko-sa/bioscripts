#!/bin/bash

# fastp trimming was superior when working with rRNAdepletion libraries
# $1 - left
# $2 - right 
# $3 - sample

#SBATCH --partition=short        # Partition (queue)
#SBATCH --time=05:00:00              # Runtime in D-HH:MM format
#SBATCH --job-name=schumacher            # Job name
#SBATCH -c 10
#SBATCH --mem=10G            # Memory needed per CPU
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL             # Type of email notification (BEGIN, END, FAIL, ALL)

fastp -i $1 -I $2 -o $3.trimmed_1.fq.gz -O $3.trimmed_2.fq.gz -h $3.fastp_html -g -x -5 -3 -w 10

