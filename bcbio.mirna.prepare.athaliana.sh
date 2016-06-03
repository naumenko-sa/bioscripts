#!/bin/bash

export APATH=/hpf/largeprojects/ccmbio/naumenko/mirna_athaliana

bcbio_setup_genome.py -f $APATH/GCF_000001735.3_TAIR10_genomic.fna \
  --gff3 -g $APATH/GCF_000001735.3_TAIR10_genomic.gff -n athaliana -b TAIR10 --mirbase ath -i star \
  --srna_gtf $APATH/ath.gtf -c 5
