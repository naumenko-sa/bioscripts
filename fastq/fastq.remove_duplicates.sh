#!/bin/bash

prinseq-lite.pl -fastq $1 -log log -derep 1 -out_good $1.derep
 