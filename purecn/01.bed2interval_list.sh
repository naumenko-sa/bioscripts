#!/bin/bash

gatk BedToIntervalList \
-I panel.bed \
-O panel.interval_list \
-SD /data/bcbio/genomes/Hsapiens/hg38/seq/hg38.dict
