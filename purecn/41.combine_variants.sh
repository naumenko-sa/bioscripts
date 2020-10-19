#!/bin/bash -l

#SBATCH --partition=priority       # Partition (queue) priority
#SBATCH --time=10:00:00            # Runtime in D-HH:MM format, 10:00:00 for hours
#SBATCH --job-name=purecn          # Job name
#SBATCH -c 1                       # cores
#SBATCH --mem=50G                  # Memory needed per CPU or --mem-per-cpu
#SBATCH --output=project_%j.out    # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err     # File to which STDERR will be written, including job ID
#SBATCH --mail-type=NONE           # Type of email notification (BEGIN, END, FAIL, ALL)


# SBATCH --job-name=purecn
# SBATCH --mem=20G
# SBATCH --export=ALL
# SBATCH -t 7-50:00
# SBATCH -p core -n 1

# $1 = normals.list of vcf.gz

# vcf requirements:
# - AD format field alt/ref
# - min genotypes 3

. .profile

date

vcf_files=""
for f in `cat $1`
do 
    vcf_files="$vcf_files -V $f"
done

gatk3 -Xmx12g \
-T CombineVariants \
--minimumN 3 \
-R $bcbio/genomes/Hsapiens/hg38/seq/hg38.fa \
-o snv.pon.vcf \
$vcf_files

bgzip snv.pon.vcf
tabix snv.pon.vcf.gz

date
