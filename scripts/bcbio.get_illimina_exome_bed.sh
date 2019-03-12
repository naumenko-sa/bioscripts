#!/bin/bash

#download and prepare exome coordinates for Illumina nextera method (bed file)
#for use in bcbio variant calling pipeline

wget http://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_nextera/nexterarapidcapture/nexterarapidcapture_exome_targetedregions_v1.2.bed
cat nexterarapidcapture_exome_targetedregions_v1.2.bed | sed s/chr//g | sed s/M/MT/g > nexterarapidcapture_exome_targetedregions_v1.2.no_chr.MT.bed