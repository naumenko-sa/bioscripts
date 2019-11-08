#!/bin/bash -l

#SBATCH --partition=short        # Partition (queue)
#SBATCH --time=10:00:00              # Runtime in D-HH:MM format
#SBATCH --job-name=bcl2fastq       # Job name
#SBATCH -c 8			    # cores
#SBATCH --mem-per-cpu=3G           # Memory needed per CPU or --mem
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL             # Type of email notification (BEGIN, END, FAIL, ALL)

date

module load bcl2fastq/2.20.0.422

bcl2fastq \
--adapter-stringency 0.9 \
--barcode-mismatches 0 \
--fastq-compression-level 4 \
--ignore-missing-bcls \
--ignore-missing-filter \
--ignore-missing-positions \
--min-log-level INFO \
--minimum-trimmed-read-length 1 \
--sample-sheet $1 \
--runfolder-dir $2 \
--output-dir $2 \
--processing-threads 8 \
--mask-short-adapter-reads 0 \
--use-bases-mask y*,y*,y*,y*

#https://github.com/indrops/indrops
#bcl2fastq --use-bases-mask y*,y*,y*,y* --mask-short-adapter-reads 0 --minimum-trimmed-read-length 0

date
