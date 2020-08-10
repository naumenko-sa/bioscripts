#!/bin/bash -l

#SBATCH --job-name=bcbio
#SBATCH --mem=50G
#SBATCH --export=ALL
#SBATCH -t 7-50:00
#SBATCH -p core -n 8

date

# run Mutect2 in T only mode to produce calls for PON
# or TO calls without matched normals

# $1 = sample_T.bam
# $2 = panel.interval_list
# $3 = snv.pon.vcf.gz

# -tumor is deprecated

. .profile

which gatk

bname=`basename $1 .bam`

gatk Mutect2 \
-R $bcbio/genomes/Hsapiens/hg38/seq/hg38.fa \
-I $1 \
-O $bname.vcf.gz \
--max-mnp-distance 0 \
--intervals $2 \
--interval-padding 50 \
--germline-resource $PURECN/af-only-gnomad.hg38.vcf.gz \
--genotype-germline-sites \
--native-pair-hmm-threads 8 \
-pon $3

tabix -f $bname.vcf.gz

date
