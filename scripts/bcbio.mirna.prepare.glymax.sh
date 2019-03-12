#!/bin/bash

bcbio_setup_genome.py -c 5 -f /hpf/largeprojects/ccmbio/naumenko/mirna/GCF_000004515.4_Glycine_max_v2.0_genomic.fna --gff3 -g \
 /hpf/largeprojects/ccmbio/naumenko/mirna/GCF_000004515.4_Glycine_max_v2.0_genomic.gff -n glycinemax -b glymax2 --mirbase gma -i star \
 --srna_gtf /hpf/largeprojects/ccmbio/naumenko/mirna/GCF_000004515.4_Glycine_max_v2.0_genomic.gtf
