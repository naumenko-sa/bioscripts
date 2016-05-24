#!/bin/bash

KPATH=/home/naumenko/work/mirna/input

bcbio_nextgen.py -w template potato.template.yaml srna12.csv \
$KPATH/HI.3550.004.RPI3.R4217_R1.fastq.gz \
$KPATH/HI.3550.004.RPI1.R4215_R1.fastq.gz \
$KPATH/HI.3550.004.RPI4.R4218_R1.fastq.gz \
$KPATH/HI.3550.004.RPI2.R4216_R1.fastq.gz \
$KPATH/HI.3550.004.RPI6.R4220_R1.fastq.gz \
$KPATH/HI.3550.004.RPI5.R4219_R1.fastq.gz \
$KPATH/HI.3550.005.RPI41.R4228_R1.fastq.gz \
$KPATH/HI.3550.005.RPI7.R4226_R1.fastq.gz \
$KPATH/HI.3550.005.RPI43.R4230_R1.fastq.gz \
$KPATH/HI.3550.005.RPI42.R4229_R1.fastq.gz \
$KPATH/HI.3550.005.RPI8.R4227_R1.fastq.gz \
$KPATH/HI.3550.005.RPI44.R4231_R1.fastq.gz 
