#!/bin/bash

date
hostname

if [ $# -lt 4 ]
then
    echo "Usage : blastn.sh qry.fasta base.fasta num_hits evalue"
    echo "Base should be ready"
    exit 0
fi

module load blast/2.6.0+
export BLASTDB=/n/data1/cores/bcbio/PIs/george_church/church2020_ecoli_high_quality_multiple_sequence_alignments_hbc04001/data/00_blastdb

qry=$1
base=$2
num_hits=$3
ev=$4

#blastn -num_threads 20 \
psiblast \
-query $qry \
-db $base -out ${qry}_vs_${base}.blastn.${ev} -evalue $ev -outfmt "6 qseqid sseqid length pident qstart qend sstart send evalue bitscore mismatch" \
-num_alignments $num_hits

date
