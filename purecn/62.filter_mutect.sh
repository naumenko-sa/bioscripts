#!/bin/bash -l

# SBATCH --partition=short       # Partition (queue) priority
# SBATCH --time=10:00:00            # Runtime in D-HH:MM format, 10:00:00 for hours
# SBATCH --job-name=filter          # Job name
# SBATCH -c 1                       # cores
# SBATCH --mem=10G                  # Memory needed per CPU or --mem-per-cpu
# SBATCH --output=project_%j.out    # File to which STDOUT will be written, including job ID
# SBATCH --error=project_%j.err     # File to which STDERR will be written, including job ID
# SBATCH --mail-type=NONE           # Type of email notification (BEGIN, END, FAIL, ALL)

#SBATCH --job-name=bcbio
#SBATCH --mem=50G
#SBATCH --export=ALL
#SBATCH -t 7-50:00
#SBATCH -p core -n 8

date

# run Mutect2 in T only mode to produce calls for PON
# or TO calls without matched normals

# $1 = unfiltered.vcf.gz

# -tumor is deprecated

. .profile

which gatk

bname=`basename $1 .vcf.gz`

gatk FilterMutectCalls \
-R $bcbio/genomes/Hsapiens/hg38/seq/hg38.fa \
-V $1 \
-O $bname.filtered.vcf.gz

tabix -f $bname.vcf.gz

date
