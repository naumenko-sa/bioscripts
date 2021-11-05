#!/bin/bash -l

#SBATCH --job-name=arriba
#SBATCH --mem=30G
#SBATCH --export=ALL
#SBATCH -t 7-50:00
#SBATCH -p core -n 1

ml bcbio-nextgen/1.2.8

#reference=/projects/ngs/reference/genomes/Hsapiens/hg19
reference=$2

date

arriba \
-x $1 \
-g ${reference}/rnaseq/ref-transcripts.gtf \
-a ${reference}/seq/hg19.fa \
-o fusions.tsv \
-O fusions.discarded.tsv \
-T -P -i chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
-b ${reference}/rnaseq/fusion-blacklist/arriba-blacklist.tsv.gz

date
