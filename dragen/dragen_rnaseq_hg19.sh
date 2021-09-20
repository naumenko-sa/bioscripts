#!/bin/bash

# build hash

#dragen \
#--build-hash-table true \ 
#--ht-build-rna-hashtable true \
#--output-directory /staging/reference/genomes/Hsapiens/hg19 \
#--ht-reference /staging/reference/genomes/Hsapiens/hg19/hg19.fa \
#--ht-alt-liftover /opt/edico/liftover/hg19_alt_liftover.sam 

date

sample=$1

mkdir dragen_${sample}

dragen \
--output-directory /staging/2021-07-29_rnaseq/dragen_${sample} \
-r /staging/reference/genomes/Hsapiens/hg38_ERCC \
--output-file-prefix $sample \
-1 /staging/2021-07-29_rnaseq/input/${sample}_R1_001.fastq.gz \
-2 /staging/2021-07-29_rnaseq/input/${sample}_R2_001.fastq.gz \
--enable-rna true \
--annotation-file /staging/reference/genomes/Hsapiens/hg38_ERCC/gencode.v38.annotation.TSL_le_3.gtf \
--enable-rna-quantification true \
--rna-quantification-library-type A \
--rna-quantification-gc-bias true \
--RGID $sample \
--RGSM $sample \
--intermediate-results-dir /staging/tmp \
--enable-rna-gene-fusion true

date
