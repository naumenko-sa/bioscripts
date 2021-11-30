#!/bin/bash

sample=M0253

mkdir /staging/$sample

dragen -f \
--tumor-fastq1 /path/${sample}_R1.fq.gz \
--tumor-fastq2 /path/${sample}_R2.fq.gz \
--enable-vcf-compression=true \
--qc-coverage-region-1=/path/M0253.bed \
--qc-coverage-filters-1 'mapq<10,bq<30' \
--enable-map-align=true \
--max-base-quality=63 \
--enable-sort=true \
--generate-sa-tags true \
--enable-bam-indexing true \
--enable-map-align-output true \
--umi-min-supporting-reads 2 \
--enable-save-bed-file true \
--umi-metrics-interval-file=/path/M0253.bed \
--umi-library-type random-simplex \
--umi-min-map-quality 1 \
--umi-max-read-error-rate 0.025 \
--umi-fuzzy-window-size 0 \
--ref-dir /staging/reference/genomes/Hsapiens/hg19 \
--output-directory /path/output \
--vc-enable-umi-liquid true \
--enable-variant-caller true \
--vc-target-bed /path/M0253.bed \
--vc-target-bed-padding 100 \
--enable-sv false \
--output-file-prefix M0253 \
--intermediate-results-dir /staging/M0253 \
--vc-enable-phasing false \
--RGID-tumor $sample \
--RGSM-tumor $sample 

# --sv-calls-regions-bed
