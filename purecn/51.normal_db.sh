#!/bin/bash -l

#SBATCH --partition=priority       # Partition (queue) priority
#SBATCH --time=10:00:00            # Runtime in D-HH:MM format, 10:00:00 for hours
#SBATCH --job-name=purecn          # Job name
#SBATCH -c 1                       # cores
#SBATCH --mem=50G                  # Memory needed per CPU or --mem-per-cpu
#SBATCH --output=project_%j.out    # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err     # File to which STDERR will be written, including job ID
#SBATCH --mail-type=NONE           # Type of email notification (BEGIN, END, FAIL, ALL)

# SBATCH --job-name=purecn
# SBATCH --mem=20G
# SBATCH --export=ALL
# SBATCH -t 7-50:00
# SBATCH -p core -n 1

date

. .profile
# 5.2 create normal.db
# $1 = snv.pon.vcf.gz
# $2 = coverage.list of coverage_loess.txt.gz gcnormalized coverage of normals
which Rscript

#Rscript $PURECN/NormalDB.R --help

Rscript $PURECN/NormalDB.R \
--outdir . \
--coveragefiles $2 \
--normal_panel $1 \
--genome hg38 \
--force

date
