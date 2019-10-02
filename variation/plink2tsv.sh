#!/bin/bash

# convert plink SNV files for a gene to tsv before importing in R
# https://www.cog-genomics.org/plink/2.0/data#export

# Options:
# $1 = dataset
# $2 = chip
# $3 = gene_name
# $4 = gene_coordinates

# Usage: 
# cd to the dir with plink results, i.e.
# cd .... dataset/data/MEGA
# 01.plink2tsv.sh [options]

# Output:
# variants: dataset.chip.gene_name.tsv
# annotation: dataset.chip.gene_name.snpinfo.tsv
# samples: dataset.chip.samples.csv

prefix=$1.$2
gene_name=$3
gene_coordinates=$4

module load plink2
plink2 --bfile data --export vcf 'id-paste=iid' --out $prefix

bgzip $prefix.vcf
tabix $prefix.vcf.gz

bcftools view -r $gene_coordinates -Oz -o $prefix.$gene_name.vcf.gz $prefix.vcf.gz

gunzip -c $prefix.$gene_name.vcf.gz | grep CHROM | sed s/"#"// | awk '{print "dataset\tchip\t"$0}'> $prefix.$gene_name.tsv
gunzip -c $prefix.$gene_name.vcf.gz | grep -v "^#" | awk -v dataset=$1 -v chip=$2 '{print dataset"\t"chip"\t"$0}' >> $prefix.$gene_name.tsv

echo "SNPname	Chromosome	Position	ReferenceAllele	AlternativeAllele	Strand	CallRateBin	RSID	Annotation" > $prefix.$gene_name.snpinfo.tsv
cat snpinfo.tsv | grep $gene_name >> $prefix.$gene_name.snpinfo.tsv

echo sample_name > $prefix.samples.csv
bcftools query -l $prefix.$gene_name.vcf.gz >> $prefix.samples.csv
