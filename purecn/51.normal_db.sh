#!/bin/bash -l

#SBATCH --job-name=purecn
#SBATCH --mem=20G
#SBATCH --export=ALL
#SBATCH -t 7-50:00
#SBATCH -p core -n 1

date

. .profile
# 5.2 create normal.db
# $1 = snv.pon.vcf.gz
# $2 = coverage.list of coverage_loess.txt.gz gcnormalized coverage of normals
which Rscript

#Rscript $PURECN/NormalDB.R --help

Rscript $PURECN/NormalDB.R \
--outdir . \
--coveragefiles $2 \
--normal_panel $1 \
--genome hg38 \
--force

date
