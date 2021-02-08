#!/bin/bash -l

# #SBATCH --partition=short          # Partition (queue) priority
# #SBATCH --time=10:00:00            # Runtime in D-HH:MM format, 10:00:00 for hours
# #SBATCH --job-name=purecn          # Job name
# #SBATCH -c 1                       # cores
# #SBATCH --mem=10G                  # Memory needed per CPU or --mem-per-cpu
# #SBATCH --output=project_%j.out    # File to which STDOUT will be written, including job ID
# #SBATCH --error=project_%j.err     # File to which STDERR will be written, including job ID
# #SBATCH --mail-type=NONE           # Type of email notification (BEGIN, END, FAIL, ALL)


#SBATCH --job-name=purecn
#SBATCH --mem=10G
#SBATCH --export=ALL
#SBATCH -t 7-50:00
#SBATCH -p core -n 1


# once OOM killed with 20G and bootstrap500
date

. .profile
which Rscript

# $1 = sample-ready_coverage_loess.txt.gz
# $2 = sample-batch-mutect2-annotated.vcf.gz

SAMPLEID=`echo $2 | awk -F '-' '{print $1}'`

echo $SAMPLEID

Rscript \
$PURECN/PureCN.R \
--out $SAMPLEID \
--tumor $1 \
--sampleid $SAMPLEID \
--vcf $2 \
--normaldb normalDB_hg38.rds \
--mappingbiasfile mapping_bias_hg38.rds \
--intervals panel.txt \
--snpblacklist $BCBIO/genomes/Hsapiens/hg38/coverage/problem_regions/repeats/simple_repeat.bed \
--genome hg38 \
--force \
--postoptimize \
--seed 123 \
--bootstrapn 500 \
--cores 1 \
--funsegmentation PSCBS \
--minpurity 0.1 --minaf 0.01 --error 0.0005

#--statsfile # mutect1

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
