#!/bin/bash

# $1 = interval_list

### 1.3 Create genomics.db
gatk GenomicsDBImport \
-R /data/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa \
-L $1 \
--genomicsdb-workspace-path snv_pon.db \
-V 12669_842__N_FFPE.chr22.vcf.gz
-V 1397_16_N_FFPE.chr22.vcf.gz
-V 2070_17_N_FFPE.chr22.vcf.gz
-V 503_17_N_FFPE.chr22.vcf.gz
-V 9855_1196__N_FFPE.chr22.vcf.gz
-V NonTNBC_22_BREAST_N_1.chr22.vcf.gz
-V SOCHG_01_OVARIAN_N.chr22.vcf.gz
-V SOCHG_02_OVARIAN_N.chr22.vcf.gz

### 1.4 Download gnomad-af only vcf:
#wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz

### 1.5. combine calls into PON.vcf

gatk CreateSomaticPanelOfNormals \
-R /data/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa \
--germline-resource af-only-gnomad.hg38.vcf.gz \
-V gendb://snv_pon_db \
-O snv_pon.vcf.gz

