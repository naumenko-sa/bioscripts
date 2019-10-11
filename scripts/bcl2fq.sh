#!/bin/bash

date

module load bcl2fastq/2.18.0.12

# Harvard Single Cell core script - no trimming of reads
# we use bcl2fastq version 2.18.0.12
bcl2fastq --use-bases-mask y*,y*,y*,y* --mask-short-adapter-reads 0 --minimum-trimmed-read-length 0

# Bauer core 
# bcl2fastq \
# --adapter-stringency 0.9 
# --barcode-mismatches 0 
# --fastq-compression-level 4 
# --ignore-missing-bcls 
# --ignore-missing-filter 
# --ignore-missing-positions 
# --min-log-level INFO 
# --minimum-trimmed-read-length 0 
# --sample-sheet /source/190531_A00794_0025_AHKKFJDMXX/SampleSheet.csv 
# --runfolder-dir /source/190531_A00794_0025_AHKKFJDMXX 
# --output-dir /output/190531_A00794_0025_AHKKFJDMXX 
# --processing-threads 8 --mask-short-adapter-reads 0 --use-bases-mask y*,y*,y*,y*

date