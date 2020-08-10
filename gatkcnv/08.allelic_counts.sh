#!/bin/bash

#SBATCH --partition=short        # Partition (queue)
#SBATCH --time=5:00:00              # Runtime in D-HH:MM format
#SBATCH --job-name=gatkcnv           # Job name
#SBATCH -c 1			    # cores
#SBATCH --mem=10G                   # Memory needed per CPU or --mem
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=NONE             # Type of email notification (BEGIN, END, FAIL, ALL)

# $1 = bam
# $2 = interval list

. .profile

bname=`basename $1 .bam`



gatk --java-options "-Xmx3g" CollectAllelicCounts \
    -L $2 \
    -I $1 \
    -R $3 \
    -O $bname.allelicCounts.tsv
