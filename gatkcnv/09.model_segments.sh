#!/bin/bash

#SBATCH --job-name=bcbio
#SBATCH --mem=10G
#SBATCH --export=ALL
#SBATCH -t 7-50:00
#SBATCH -p core -n 1

# $1 = T.denoisedCR.tsv

bname=`echo $1 | awk -F "." '{print $1}'`

gatk --java-options "-Xmx4g" ModelSegments \
--denoised-copy-ratios $1 \
--output . \
--output-prefix $bname
