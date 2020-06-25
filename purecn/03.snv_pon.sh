#!/bin/bash

ls -1 *.for_pon.vcf.gz > sample_list.txt

gatk CreateSomaticPanelOfNormals \
-vcfs sample_list.txt \
-O snv_pon.vcf.gz
