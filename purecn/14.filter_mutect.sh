#!/bin/bash -l

#SBATCH --job-name=bcbio
#SBATCH --mem=50G
#SBATCH --export=ALL
#SBATCH -t 7-50:00
#SBATCH -p core -n 8

date

# run Mutect2 in T only mode to produce calls for PON
# or TO calls without matched normals

# $1 = unfiltered.vcf.gzsample_N.bam

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
