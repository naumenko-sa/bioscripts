#!/bin/bash

# http://lindenb.github.io/jvarkit/BamStats05.html
# bamstats05 - groups coverage by gene
# really fast: 2min for a big bam file and 100 genes
# arguments: bam and bed
# by default uses protein coding genes

# pull exone coordinates with genes.R and do bedtools sort | bedtools merge -c 4 -o first 

# https://slurm.schedmd.com/sbatch.html
# https://wiki.rc.hms.harvard.edu/display/O2

#SBATCH --partition=short        # Partition (queue) priority
#SBATCH --time=1:00:00              # Runtime in D-HH:MM format, 10:00:00 for hours
#SBATCH --job-name=project          # Job name
#SBATCH -c 5			    # cores
#SBATCH --mem=15G           # Memory needed per CPU or --mem-per-cpu
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=             # Type of email notification (BEGIN, END, FAIL, ALL)

JVARKIT_PATH=~/Desktop/tools/JVARKIT/jvarkit.jar

# if you need all reads, add -f "" - empty filter, by default it filters out duplicated reads

# could be gz'ed
input_file=$1
bed=$2
ref=$3
mapq=$4

if [ -z $bed ]
then
    bed=~/cre/data/protein_coding_genes.exons.fixed.sorted.bed
fi

if [[ ! -n "$mapq" ]]; then
    mapq=30
fi

ext="${input_file#*.}"
sample="${input_file%%.*}"

# note stringent mapq and mincov
if [[ "$ext" == "bam" ]];then

    /usr/bin/java -Xmx10G -jar $JVARKIT_PATH bamstats05 \
        --bed $bed \
        --mapq $mapq \
        --mincoverage 10 $input_file > $sample.mq30.minc10.coverage.tsv

elif [[ "$ext" == "cram" ]];then

    java -Xmx20G -jar $JVARKIT_PATH bamstats05 \
        --bed $bed \
        --mapq $mapq \
        --mincoverage 10 \
        --reference $ref $input_file > $sample.mq${mapq}.minc10.coverage.tsv
fi

