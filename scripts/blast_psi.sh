#!/bin/bash

#SBATCH --partition=short           # Partition (queue) priority
#SBATCH --time=12:00:00             # Runtime in D-HH:MM format, 10:00:00 for hours
#SBATCH --job-name=blast            # Job name
#SBATCH -c 20			    # cores
#SBATCH --mem=20G                  # total Memory or use --mem-per-cpu
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL             # Type of email notification (BEGIN, END, FAIL, ALL)

date
echo $1
echo $2

if [ $# -lt 4 ]
then
    echo "Usage : blastn.sh qry.fasta base.fasta num_hits evalue"
    echo "Base should be ready"
    exit 0
fi

# set BLASTDB variable
# only can do negative_taxidlist, not both lists

qry=$1
base=$2
num_hits=$3
ev=$4

psiblast -num_threads 20 \
-query $qry \
-db $base \
-out ${qry}_vs_${base}.psiblast.num_hits_${num_hits}.ev_${ev} \
-evalue $ev \
-outfmt "6 qseqid sseqid length qlen slen pident qstart qend sstart send evalue bitscore mismatch gaps staxids sscinames scomnames" \
-num_alignments $num_hits

date
