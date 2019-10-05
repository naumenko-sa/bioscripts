#!/bin/bash -l

#SBATCH --job-name=bcbio
#SBATCH --mem=20G
#SBATCH --export=ALL
#SBATCH -t 10:00:00
#SBATCH -p core -n 1

date

# parse duplex barcodes

module use /projects/ngs/local/software/modules
module load bcbio-nextgen/latest-devel

which bcbio_fastq_umi_prep.py
bcbio_fastq_umi_prep.py autopair --tag1 8 --tag2 8 $1 $2

date
