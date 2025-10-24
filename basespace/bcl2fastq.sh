#!/bin/bash

# install via conda
# https://anaconda.org/dranew/bcl2fastq

# NovaSeqX requires bcl-convert
# rpm2cpio bcl-convert-4.1.7-2.el7.x86_64.rpm | cpio -idmv

# needs sudo for some reason to write files
# input dir is a run dir
# to alter sample sheet it has to be copied
# sudo /data/tools/bcl-convert --bcl-input-directory $1 --output-directory ./fastq --sample-sheet ./SampleSheet.csv &> log.txt
