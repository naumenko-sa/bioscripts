#!/bin/bash

export PYTHONPATH=
export LD_LIBRARY_PATH=

/usr/bin/python bcbio_nextgen_install.py /hpf/largeprojects/ccmbio/naumenko/tools/bcbio --tooldir=/hpf/largeprojects/ccmbio/tools/bcbio --genomes GRCh37 --aligners bwa --datatarget vep
