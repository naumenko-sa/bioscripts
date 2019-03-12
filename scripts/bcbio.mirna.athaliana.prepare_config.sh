#!/bin/bash

KPATH=/hpf/largeprojects/ccmbio/naumenko/mirna2/input

bcbio_nextgen.py -w template bcbio.mirna.athaliana.template.yaml bcbio.mirna.athaliana.35samples.csv \
$KPATH/12T1R1.MPS12300453-E06.3550.3.single.fastq.gz \
$KPATH/12T1R2.MPS12300453-F06.3550.2.single.fastq.gz \
$KPATH/12T1R3.MPS12302114-F09.3550.3.single.fastq.gz \
$KPATH/12T2R1.MPS12302118-D10.3550.3.single.fastq.gz \
$KPATH/12T2R2.MPS12300453-H06.3550.3.single.fastq.gz \
$KPATH/12T2R3.MPS12302118-F10.3550.3.single.fastq.gz \
$KPATH/12T3R1.MPS12302118-G10.3550.3.single.fastq.gz \
$KPATH/12T3R2.MPS12302118-H10.3550.3.single.fastq.gz \
$KPATH/12T3R3.MPS12300333-A10.3550.2.single.fastq.gz \
$KPATH/12T4R1.MPS12300333-B10.3550.3.single.fastq.gz \
$KPATH/12T4R2.MPS12300333-C10.3550.2.single.fastq.gz \
$KPATH/12T4R3.MPS12300333-D10.3550.2.single.fastq.gz \
$KPATH/1T1R1.MPS12302118-A07.3550.2.single.fastq.gz \
$KPATH/1T1R2.MPS12300453-A06.3550.3.single.fastq.gz \
$KPATH/1T1R3.MPS12302118-C07.3550.2.single.fastq.gz \
$KPATH/1T2R1.MPS12302118-D07.3550.2.single.fastq.gz \
$KPATH/1T2R2.MPS12302118-E07.3550.2.single.fastq.gz \
$KPATH/1T2R3.MPS12302114-C09.3550.2.single.fastq.gz \
$KPATH/1T3R1.MPS12302118-G07.3550.2.single.fastq.gz \
$KPATH/1T3R2.MPS12302118-H07.3550.2.single.fastq.gz \
$KPATH/1T3R3.MPS12300333-A09.3550.3.single.fastq.gz \
$KPATH/1T4R2.MPS12300333-C09.3550.3.single.fastq.gz \
$KPATH/1T4R3.MPS12300333-D09.3550.3.single.fastq.gz \
$KPATH/6T1R1.MPS12302118-E08.3550.2.single.fastq.gz \
$KPATH/6T1R2.MPS12302118-F08.3550.2.single.fastq.gz \
$KPATH/6T1R3.MPS12302118-G08.3550.2.single.fastq.gz \
$KPATH/6T2R1.MPS12302118-H08.3550.2.single.fastq.gz \
$KPATH/6T2R2.MPS12302118-A09.3550.2.single.fastq.gz \
$KPATH/6T2R3.MPS12302118-B09.3550.2.single.fastq.gz \
$KPATH/6T3R1.MPS12302118-C09.3550.2.single.fastq.gz \
$KPATH/6T3R2.MPS12302118-D09.3550.3.single.fastq.gz \
$KPATH/6T3R3.MPS12300333-E09.3550.3.single.fastq.gz \
$KPATH/6T4R1.MPS12300333-F09.3550.3.single.fastq.gz \
$KPATH/6T4R2.MPS12302114-E09.3550.3.single.fastq.gz \
$KPATH/6T4R3.MPS12300333-H09.3550.3.single.fastq.gz
