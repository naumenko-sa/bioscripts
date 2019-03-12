#!/bin/bash
date

if [ $# -lt 2 ];then
    echo "Usage : goblastp.sh qry.fasta base.fasta num_hits evalue";
    echo "Base should be ready";
    exit 0;
fi;

qry=$1
base=$2
num_hits=$3
ev=$4

/mnt/lustre/tools/ncbi-blast-2.2.29+/bin/blastp -task blastp -num_threads 24 -query $qry -db $base -out ${qry}_vs_${base}.blastp.${ev}.xml -evalue $ev -outfmt 5 -num_alignments $num_hits

date

