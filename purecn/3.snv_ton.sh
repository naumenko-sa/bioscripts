#!/bin/bash -l

#SBATCH --job-name=bcbio
#SBATCH --mem=50G
#SBATCH --export=ALL
#SBATCH -t 7-50:00
#SBATCH -p core -n 8

date

module load bcbio-nextgen/latest-devel

# $1 = sample_T.bam
# $2 = panel.interval_list
# $3 = snv.pon.vcf.gz

# -tumor is deprecated

which gatk

bcbio=/projects/ngs/reference/UpdateGenomesBcbio
germline_resource=/projects/ngs/local/users/kmhr378/2020-06-30_DEV1534_purecn/02_reference/af-only-gnomad.hg38.vcf.gz

bname=`basename $1 .bam`

gatk Mutect2 \
-R /$bcbio/Hsapiens/hg38/seq/hg38.fa \
-I $1 \
-O $bname.vcf.gz \
--max-mnp-distance 0 \
--intervals $2 \
--interval-padding 50 \
--germline-resource $germline_resource \
--genotype-germline-sites \
--native-pair-hmm-threads 16 \
--panel-of-normals $3

tabix -f $bname.vcf.gz

date
