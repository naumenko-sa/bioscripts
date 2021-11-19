#!/bin/bash

# https://slurm.schedmd.com/sbatch.html
# https://wiki.rc.hms.harvard.edu/display/O2

#SBATCH --partition=short           # Partition (queue) priority
#SBATCH --time=12:00:00             # Runtime in D-HH:MM format, 10:00:00 for hours
#SBATCH --job-name=blast            # Job name
#SBATCH -c 20			    # cores
#SBATCH --mem=50G                    # total Memory or use --mem-per-cpu
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=NONE            # Type of email notification (BEGIN, END, FAIL, ALL)

# makeblastdb -in pol.fasta -dbtype 'nucl'
# https://www.ncbi.nlm.nih.gov/nuccore/JX198552.1?report=fasta
# ~ 5 min for 200k sequences
date

echo blast_xml

ml blast/2.10.1+

qry=$1
base=$2

# default settings are suggested
# 1
# num_hits=$3
# 0.01
# ev=$4

rm ${qry}_vs_${base}.xml

blastn -num_threads 20 \
-query $qry \
-db $base \
-out ${qry}_vs_${base}.xml \
-dust no \
-num_alignments 1 \
-outfmt 5
#-evalue $ev \
#-num_alignments $num_hits

#-outfmt "6 qseqid sseqid length qlen slen pident qstart qend sstart send evalue bitscore mismatch gaps staxids sscinames scomnames" \
#-num_alignments $num_hits

date
