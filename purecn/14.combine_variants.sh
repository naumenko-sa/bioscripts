#!/bin/bash -l

#SBATCH --job-name=purecn
#SBATCH --mem=29G
#SBATCH --export=ALL
#SBATCH -t 7-50:00
#SBATCH -p core -n 1

date
bcbio=/projects/ngs/reference/UpdateGenomesBcbio

vcf_files=""
for f in *.vcf.gz
do 
    vcf_files="$vcf_files -V $f"
done

gatk3 -Xmx12g \
-T CombineVariants \
--minimumN 3 \
-R $bcbio/Hsapiens/hg38/seq/hg38.fa \
-o snv_pon.vcf \
$vcf_files

bgzip snv_pon.vcf
tabix snv_pon.vcf.gz
date