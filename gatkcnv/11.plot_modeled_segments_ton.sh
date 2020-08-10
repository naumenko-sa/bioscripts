#!/bin/bash

#SBATCH --partition=short        # Partition (queue)
#SBATCH --time=5:00:00              # Runtime in D-HH:MM format
#SBATCH --job-name=gatkcnv           # Job name
#SBATCH -c 1			    # cores
#SBATCH --mem=10G                   # Memory needed per CPU or --mem
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=NONE             # Type of email notification (BEGIN, END, FAIL, ALL)

# $1 = T.denoisedCR.tsv
# $2 hg38.dict

. .profile

bname=`basename $1 .denoisedCR.tsv`

gatk PlotModeledSegments \
--denoised-copy-ratios $1 \
--segments $bname.modelFinal.seg \
--sequence-dictionary $bcbio/genomes/Hsapiens/hg38/seq/hg38.dict \
--minimum-contig-length 46709983 \
--output segment_plots \
--output-prefix $bname
