#!/bin/bash
date
hostname

if [ $# -lt 4 ]
then
    echo "Usage : blastn.sh qry.fasta base.fasta num_hits evalue"
    echo "Base should be ready"
    exit 0
fi

qry=$1
base=$2
num_hits=$3
ev=$4

module load ncbi

blastn -task blastn -query $qry -db $base -out ${qry}_vs_${base}.blastn.${ev} -evalue $ev -outfmt "6 qseqid sseqid length pident qstart qend sstart send evalue bitscore mismatch" -num_alignments $num_hits -num_threads 20

date

