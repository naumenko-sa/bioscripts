#!/bin/bash

# http://lindenb.github.io/jvarkit/BamStats05.html
# uses bamstats05 - groups coverage by gene
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
#SBATCH --mail-type=ALL             # Type of email notification (BEGIN, END, FAIL, ALL)

# wants java1.8

BAMSTATS_PATH=~/work/tools/jvarkit/dist

# if you need all reads, add -f "" - empty filter, by default it filters out some duplicated reads

bed=$2
bam=$1

if [ -z $bed ]
then
    bed=~/cre/data/protein_coding_genes.exons.fixed.sorted.bed
fi

java -Xmx10G -jar ${BAMSTATS_PATH}/bamstats05.jar \
    -B $bed \
    $bam > $bam.coverage
