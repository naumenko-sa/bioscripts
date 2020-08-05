#!/bin/bash

# in createsomaticpanelofnormals
# -V does not work for multiple samples
# -vcfs --vcfs do not work
# only works through genomicsdb

ls -1 *.for_pon.vcf.gz | awk -F "." '{print $1"\t"$0}' > sample_list.tsv

gatk GenomicsDBImport \
-R /data/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa \
--intervals $1 \
--sample-name-map sample_list.tsv \
--genomicsdb-workspace-path pon_db
