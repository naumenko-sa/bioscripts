#!/bin/bash

#SBATCH --partition=short        # Partition (queue)
#SBATCH --time=5:00:00              # Runtime in D-HH:MM format
#SBATCH --job-name=gatkcnv           # Job name
#SBATCH -c 1			    # cores
#SBATCH --mem=10G                   # Memory needed per CPU or --mem
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=NONE             # Type of email notification (BEGIN, END, FAIL, ALL)

# $1 = sort.bam
# $2 = interval.list - not gc annotated
# $3 = TSV for tsv output, default is hdf5

date

bname=`basename $1 .bam`

if [ -z $3 ]
then
    gatk --java-options '-Xms500m -Xmx131859m -XX:+UseSerialGC -Djava.io.tmpdir=.' \
    CollectReadCounts \
    -I $1 \
    -L $2 \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O $bname.counts.hdf5
else
   gatk --java-options '-Xms500m -Xmx131859m -XX:+UseSerialGC -Djava.io.tmpdir=.' \
    CollectReadCounts \
    -I $1 \
    -L $2 \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O $bname.counts.tsv \
    --format TSV
fi

date
