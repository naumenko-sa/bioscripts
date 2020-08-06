#!/bin/bash -l

#SBATCH --job-name=purecn
#SBATCH --mem=20G
#SBATCH --export=ALL
#SBATCH -t 7-50:00
#SBATCH -p core -n 1

date

. /projects/ngs/local/users/kmhr378/2020-07-09_BS_benchmark/.bash_profile
bcbio=/projects/ngs/local/users/kmhr378/2020-07-09_BS_benchmark/bcbio
PURECN=$bcbio/anaconda/envs/r36/lib/R/library/PureCN/extdata

export PATH=$bcbio/anaconda/envs/r36/bin:$PATH

which Rscript

# $1 = tumor.coverage_loess.txt.gz
# $2 = mutect.vcf.gz

SAMPLEID=`echo $1 | awk -F '.' '{print $1}'`

Rscript \
$PURECN/PureCN.R \
--out . \
--tumor $1 \
--sampleid $SAMPLEID \
--vcf $2 \
--statsfile $2.stats \
--normaldb normalDB_hg38.rds \
--mappingbiasfile mapping_bias_hg38.rds \
--intervals panel.txt \
--snpblacklist hg19_simpleRepeats.bed \
--genome hg38 \
--force --postoptimize --seed 123

#Rscript $PURECN/PureCN.R \
#--sampleid $1 \
#--out $1.out \
#--tumor $1.counts.hdf5 \
#--logratiofile $1.denoisedCR.tsv \
#--segfile $1.modelFinal.seg \
#--mappingbiasfile $2 \
#--vcf $1.vcf.gz \
#--statsfile $1.vcf.gz.stats \
#--genome hg38 \
#--funsegmentation Hclust \
#--force \
#--postoptimize \
#--seed 123 \
#--genome hg38

date
