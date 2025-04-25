#!/bin/bash

# download, filter, and concatenate high coverage 3202 sample 1KG dataset from NYGenomeCentre
# no genotypes, AF, AN only
# output: 20G 3202.vcf.gz

for f in {`seq 1 22`,X,Y}
do
    if [ ! -f chr${f}.vcf.gz ];then
        wget -c https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${f}.recalibrated_variants.annotated.vcf.gz -O chr${f}.vcf.gz
        wget -c https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${f}.recalibrated_variants.annotated.vcf.gz.tbi -O chr${f}.vcf.gz.tbi
        bcftools view -f "PASS,." chr${f}.vcf.gz -Ov | bgzip -c > chr${f}.pass.vcf.gz
        tabix chr${f}.pass.vcf.gz
        rm chr${f}.vcf.gz
        rm chr${f}.vcf.gz.tbi
        mv chr${f}.pass.vcf.gz chr${f}.vcf.gz
        mv chr${f}.pass.vcf.gz.tbi chr${f}.vcf.gz.tbi
    fi
done

bcftools concat chr1.vcf.gz chr2.vcf.gz chr3.vcf.gz chr4.vcf.gz chr5.vcf.gz \
                chr6.vcf.gz chr7.vcf.gz chr8.vcf.gz chr9.vcf.gz chr10.vcf.gz \
                chr11.vcf.gz chr12.vcf.gz chr13.vcf.gz chr14.vcf.gz chr15.vcf.gz \
                chr16.vcf.gz chr17.vcf.gz chr18.vcf.gz chr19.vcf.gz chr20.vcf.gz \
                chr21.vcf.gz chr22.vcf.gz chrX.vcf.gz chrY.vcf.gz -Ov > 3202.vcf

bgzip 3202.vcf
tabix 3202.vcf.gz
