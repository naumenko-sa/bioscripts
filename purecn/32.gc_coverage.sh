#!/bin/bash -l

# SBATCH --partition=short          # Partition (queue) priority
# SBATCH --time=10:00:00            # Runtime in D-HH:MM format, 10:00:00 for hours
# SBATCH --job-name=purecn         # Job name
# SBATCH -c 8                       # cores
# SBATCH --mem=50G                  # Memory needed per CPU or --mem-per-cpu
# SBATCH --output=project_%j.out    # File to which STDOUT will be written, including job ID
# SBATCH --error=project_%j.err     # File to which STDERR will be written, including job ID
# SBATCH --mail-type=NONE           # Type of email notification (BEGIN, END, FAIL, ALL)


#SBATCH --job-name=purecn
#SBATCH --mem=20G
#SBATCH --export=ALL
#SBATCH -t 7-50:00
#SBATCH -p core -n 1

date

. .profile
### 3.2 gc_coverage
# $1 = sample.bam
# $2 = intervals.txt = output of process_intervals

# output:
# - coverage.txt.gz
# - coverages_loess.png
# - coverage_loess.txt.gz
# - coverage_loess_qc.txt

which Rscript

Rscript \
$PURECN/Coverage.R \
--outdir . \
--bam $1 \
--intervals $2

date
