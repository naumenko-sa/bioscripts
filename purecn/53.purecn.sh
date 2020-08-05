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
# $2 = list of hgf5 files gcnormalized coverage of normals
PURECN=$bcbio/anaconda/envs/r36/lib/R/library/PureCN/extdata

export PATH=$bcbio/anaconda/envs/r36/bin:$PATH

which Rscript

# $1 = sample_id
# $2 = mapping_bias.rds

Rscript $PURECN/PureCN.R \
--sampleid $1 \
--out $1.out \
--tumor ${1}-target-coverage.hdf5 \
--logratiofile $1.denoisedCR.tsv \
--segfile $1.modelFinal.seg \
--mappingbiasfile $2 \
--vcf $1.vcf.gz \
--statsfile $1.vcf.gz.stats \
--genome hg38 \
--funsegmentation Hclust \
--force \
--postoptimize \
--seed 123

date
