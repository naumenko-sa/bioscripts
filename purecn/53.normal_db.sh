#!/bin/bash -l

#SBATCH --job-name=purecn
#SBATCH --mem=20G
#SBATCH --export=ALL
#SBATCH -t 7-50:00
#SBATCH -p core -n 1

date

. /projects/ngs/local/users/kmhr378/2020-07-09_BS_benchmark/.bash_profile
bcbio=/projects/ngs/local/users/kmhr378/2020-07-09_BS_benchmark/bcbio

# 5.2 create normal.db
# $1 = snv.pon.vcf.gz
# $2 = coverage.list of coverage_loess.txt.gz gcnormalized coverage of normals
PURECN=$bcbio/anaconda/envs/r36/lib/R/library/PureCN/extdata

export PATH=$bcbio/anaconda/envs/r36/bin:$PATH

which Rscript

#Rscript $PURECN/NormalDB.R --help

Rscript $PURECN/NormalDB.R \
--outdir . \
--coveragefiles $2 \
--normal_panel $1 \
--genome hg38 \
--force

date
