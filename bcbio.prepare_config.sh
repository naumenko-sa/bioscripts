#!/bin/bash

KPATH=/home/naumenko/work/project_muscular/input

bcbio_nextgen.py -w template ~/bioscripts/bcbio.exome.template.yaml md.csv \
$KPATH/blood1_1.fq.gz  \
$KPATH/blood1_2.fq.gz  \
$KPATH/blood2_1.fq.gz  \
$KPATH/blood2_2.fq.gz  \
$KPATH/fibroblast5_1.fq.gz  \
$KPATH/fibroblast5_2.fq.gz  \
$KPATH/muscle1_1.fq.gz  \
$KPATH/muscle1_2.fq.gz  \
$KPATH/muscle2_1.fq.gz  \
$KPATH/muscle2_2.fq.gz  \
$KPATH/myotubes5_1.fq.gz  \
$KPATH/myotubes5_2.fq.gz
