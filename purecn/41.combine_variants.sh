#!/bin/bash -l

#SBATCH --job-name=purecn
#SBATCH --mem=20G
#SBATCH --export=ALL
#SBATCH -t 7-50:00
#SBATCH -p core -n 1

# $1 = normals.list of vcf.gz

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
