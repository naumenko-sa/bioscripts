#!/bin/bash

# extract reads from the table
for f in `cat samples_no_prefix.txt`; do cat all_discordant.tsv | grep "^$f" | awk -F '\t' '{if ($11 == 60 && $28 >=35 && $32 >= 20 && $25 <=1 ) {print $0}}' | grep -v chrUn | grep "On Target:Genomic" > $f.ontarget_genomic.tsv;done;
wc -l *ontarget_genomic.tsv

# read names
for f in `cat samples_no_prefix.txt`;do cat $f.ontarget_genomic.tsv | awk -F '\t' '{print $3}'> $f.ontarget_genomic.names.txt;done;

# fasta
for f in `cat samples_no_prefix.txt`;do python3  extract_reads.py --bam s_${f}-ready.bam --names $f.ontarget_genomic.names.txt --out $f.ontarget_genomic.fasta;done;
 
# parse circle seq excel to fasta
cat CIRCLE-seq\ gMH_full.csv | sed 1d | awk -F "\t" '{if (length($3)>0){print ">"$1; print$3}}' > circle_seq.fasta

# make blastb
ml BLAST+
makeblastdb -in mm10_gene.fa -dbtype 'nucl'

#blast
for f in `cat samples_no_prefix.txt`;do blastn.sh $f.ontarget_genomic.fasta mm10_gene.fa 5 0.01;done;

# proximal
for f in `cat samples_no_prefix.txt`;do python3 find_proximal.py $f.ontarget_genomic.fasta_vs_mm10_gene.fa.blastn.0.01 circle_seq.fasta_vs_mm10_gene.fa.blastn.0.01 10000 > $f.proximal.circle_seq.csv;done;

# filter_proximal
for f in `cat samples_no_prefix.txt`;do cat $f.proximal.circle_seq.csv  | grep -v chr4 | grep -v gene > $f.proximal.circle_seq.no_target.csv;done;

# statistics
for f in *no_target.csv; do echo $f `cat $f | sed 1d | wc -l`;done;

# cas offinder
cat cas_offinder.csv  | sed 1d | awk '{print ">"$2"_"$3"_"$4;print toupper($2)}' > cas_offinder.fasta
