#!/bin/bash -l

#SBATCH --partition=priority       # Partition (queue) priority
#SBATCH --time=10:00:00            # Runtime in D-HH:MM format, 10:00:00 for hours
#SBATCH --job-name=purecn          # Job name
#SBATCH -c 1                       # cores
#SBATCH --mem=20G                  # Memory needed per CPU or --mem-per-cpu
#SBATCH --output=project_%j.out    # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err     # File to which STDERR will be written, including job ID
#SBATCH --mail-type=NONE           # Type of email notification (BEGIN, END, FAIL, ALL)

## SBATCH --job-name=purecn
## SBATCH --mem=20G
## SBATCH --export=ALL
## SBATCH -t 7-50:00
## SBATCH -p core -n 1

date
. .profile

# input:
# - $1 = panel.bed
# output:
# - panel.txt

which Rscript

bname=`basename $1 .bed`

Rscript \
$PURECN/IntervalFile.R \
--infile $1 \
--fasta $bcbio/genomes/Hsapiens/hg38/seq/hg38.fa \
--outfile $bname.txt \
--offtarget \
--genome hg38 \
--export panel.optimized.bed \
--mappability $PURECN/GCA_000001405.15_GRCh38_no_alt_analysis_set_100.bw

date
