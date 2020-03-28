#!/bin/bash

#SBATCH --job-name=bcbio
#SBATCH --mem=30G
#SBATCH --export=ALL
#SBATCH -t 7-50:00
#SBATCH -p core -n 1

date

module load bcbio-nextgen/1.2.0

fgbio -Xmx1g --version

sample=$1

unset JAVA_HOME && \
fgbio \
-Xms750m -Xmx30g -XX:+UseSerialGC \
--tmp-dir . \
--async-io=true \
--compression=0 \
GroupReadsByUmi \
--edits=1 \
--min-map-q=1 \
-t RX \
-s adjacency \
-i ${sample}-umi.bam \
-f group_reads_family_histogram.tsv | \
fgbio \
-Xms750m \
-Xmx30g \
-XX:+UseSerialGC \
--tmp-dir . \
--async-io=true \
--compression=0 \
CallMolecularConsensusReads \
--min-input-base-quality=2 \
--min-reads=3 \
--max-reads=100000 \
--output-per-base-tags=false \
--sort-order=:none: \
-i /dev/stdin \
-o /dev/stdout | \
fgbio \
-Xms750m \
-Xmx30g \
-XX:+UseSerialGC \
--tmp-dir . \
--async-io=true \
--compression=0 \
FilterConsensusReads \
--min-reads=3 \
--min-base-quality=13 \
--max-base-error-rate=0.1 \
-r /projects/ngs/reference/genomes/Hsapiens/hg38/seq/hg38.fa \
-i /dev/stdin -o /dev/stdout | \
bamtofastq \
collate=1 \
T=cumi-1-bamtofastq-tmp \
F=${sample}-sort-cumi-1.fq.gz \
F2=${sample}-sort-cumi-2.fq.gz tags=cD,cM,cE gz=1
