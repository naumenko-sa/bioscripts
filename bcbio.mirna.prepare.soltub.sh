#!/bin/bash

bcbio_setup_genome.py -f /hpf/largeprojects/ccmbio/naumenko/mirna/GCF_000226075.1_SolTub_3.0_genomic.fna \
  --gff3 -g /hpf/largeprojects/ccmbio/naumenko/mirna/GCF_000226075.1_SolTub_3.0_genomic.gff -n stuberosum -b soltub3 --mirbase stu -i star \
  --srna_gtf /hpf/largeprojects/ccmbio/naumenko/mirna/GCF_000226075.1_SolTub_3.0_genomic.gtf -c 5