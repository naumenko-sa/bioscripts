#!/bin/bash

# install size ~20G
sudo docker pull ensemblorg/ensembl-vep
# install cache and reference fasta
sudo docker run -t -i -v /data/vep_data:/opt/vep/.vep:Z ensemblorg/ensembl-vep INSTALL.pl
