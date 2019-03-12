#!/bin/bash
date

if [ $# -lt 2 ];then
    echo "Usage : goblastx.sh qry.fasta base.fasta num_hits evalue";
    echo "Base should be ready";
    exit 0;
fi;

qry=$1
base=$2
num_hits=$3
ev=$4

/mnt/lustre/tools/ncbi-blast-2.2.29+/bin/blastx -num_threads 24 -query $qry -db $base -out ${qry}_vs_${base}.blastx.${ev} -evalue $ev -outfmt "6 qseqid sseqid length pident qstart qend sstart send evalue" -num_alignments $num_hits

date

