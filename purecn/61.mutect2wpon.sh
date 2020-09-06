#!/bin/bash -l

# SBATCH --partition=short       # Partition (queue) priority
# SBATCH --time=10:00:00            # Runtime in D-HH:MM format, 10:00:00 for hours
# SBATCH --job-name=mutect2pon          # Job name
# SBATCH -c 8                       # cores
# SBATCH --mem=50G                  # Memory needed per CPU or --mem-per-cpu
# SBATCH --output=project_%j.out    # File to which STDOUT will be written, including job ID
# SBATCH --error=project_%j.err     # File to which STDERR will be written, including job ID
# SBATCH --mail-type=NONE           # Type of email notification (BEGIN, END, FAIL, ALL)


#SBATCH --job-name=mutect2
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
--germline-resource $PURECN/af_only_gnomad.vcf.gz \
--genotype-germline-sites \
--native-pair-hmm-threads 8 \
-pon $3

tabix -f $bname.vcf.gz

date
