#!/bin/bash

#SBATCH --job-name=bcbio
#SBATCH --mem=10G
#SBATCH --export=ALL
#SBATCH -t 7-50:00
#SBATCH -p core -n 1

# $1 = input bam.counts.hdf5
# $2 = pon.hdf5
# $3 = panel.gcannotated.tsv - necessary for PureCN

bname=`basename $1 .counts.hdf5`

gatk --java-options "-Xmx12g" \
DenoiseReadCounts \
-I $1 \
--count-panel-of-normals $2 \
--standardized-copy-ratios $bname.standardizedCR.tsv \
--denoised-copy-ratios $bname.denoisedCR.tsv \
--annotated-intervals $3
