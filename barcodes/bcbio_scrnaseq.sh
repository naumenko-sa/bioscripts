#!/bin/bash

# step1: parse barcode information from reads 2,3,4 to the fastq read name,
# CELL_[value]:UMI_[value]:sample_[value]. umis supports many protocols, however,
# the downside is speed - this step can take up to 3-4 days::

#python umis fastqtransform \
#--separate_cb umis/harvard-indrop-v3-transform.json --cores 16  \
#project_1.fq.gz project_2.fq.gz project_3.fq.gz project_4.fq.gz | \
#seqtk seq -L 20 - | gzip > project.umitransformed.fq.gz

# step2: create a fastq file for each sample::
#python umis demultiplex_samples \
#--nedit 1 --barcodes sample_barcodes.csv \
#--out_dir demultiplexed project.umitransformed.fq.gz

# step3*::
#python umis cb_filter \
#--cores 16 --bc1 harvard-indrop-v3-cb1.txt.gz --nedit 1 \
#--bc2 harvard-indrop-v3-cb2.txt.gz demultiplexed/[sample-barcodeAATTTTT].fq  | \
#gzip -c > project-sample_barcode.filtered.fq.gz

# step 4*: create cellular barcode histrogram::
#python umis cb_histogram project-sample_barcode.filtered.fq.gz > cb-histogram.txt

# step 5: create index genome for rapmap::
#rapmap quasiindex -k 31 -i mm10 -t mm10/ref-transcripts.fa

# step 6*. align reads with rapmap::
#rapmap quasimap -t 16 -i mm10 \
#-r <(gzip -cd project-sample_barcode.filtered.fq.gz) | \
#samtools sort -@ 16 -m 1G  -T project-sample_barcode-sorttmp \
#-o project-sample_barcode.bam /dev/stdin
#samtools index -@ 16 project-sample_barcode.bam project-sample_barcode.bam.bai

# step 7*: count transcripts::
#python umis fasttagcount --cb_cutoff 1000 \
--genemap ref-transcripts-tx2gene.tsv
--cb_histogram project-sample_barcode/cb-histogram.txt \
--umi_matrix project-sample_barcode-dupes.mtx.full \
project-sample_barcode.bam project-sample_barcode.mtx.full

python umis sparse project-sample_barcode.mtx.full \
project-sample_barcode.mtx

python umis sparse project-sample_barcode-dupes.mtx.full \
project-sample_barcode-dupes.mtx

# step 8. Concatenate all cb-histogram-filtered.txt files::
cat project-[all-barcodes]/cb-histogram-filtered.txt > cb-histogram.txt
