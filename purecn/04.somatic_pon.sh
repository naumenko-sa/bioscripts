#!/bin/bash

gatk CreateSomaticPanelOfNormals \
-R /data/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa \
--germline-resource af-only-gnomad.vcf.gz \
-V gendb://pon_db \
-O snv_pon.vcf.gz
