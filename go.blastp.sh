#!/bin/bash
#blastp for ortholog pipeline

date

if [ $# -lt 2 ];then
    echo "Usage : goblastp.sh qry.fasta base.fasta";
    echo "Base should be ready";
    exit 0;
fi;

qry=$1
base=$2
ev=1e-5

/mnt/lustre/tools/ncbi-blast-2.2.29+/bin/blastp -task blastp -num_threads 23 -query $qry -db $base -out ${qry}_vs_${base}.blastp.${ev} -evalue $ev -outfmt 6
# qseqid sseqid qlen slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore

date

