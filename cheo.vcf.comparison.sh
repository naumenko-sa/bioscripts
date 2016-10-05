#!/bin/bash

#compares original vcf and the new one
#$1 - base name for files base.orig.vcf and base.bcbio.vcf

module load vcftools
module load bedtools

function prepare
{
    cat $1.bcbio.vcf | sed s/"^MT"/"^M"/ > $1.bcbio.M.vcf
  
    for f in {orig,bcbio.M};
    do 
	vcftools --vcf $1.$f.vcf --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 \
         --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 --chr X --chr Y --chr M --remove-indels --recode --recode-INFO-all --out $1.$f.SNP
    done;

    for f in {orig,bcbio.M};
    do 
	vcftools --vcf $1.$f.SNP.recode.vcf --min-meanDP 10 --recode --recode-INFO-all --out $1.$f.SNP.dp10
    done;
    
}

function lost_variants
{
    for f in {SNP,SNP.dp10};
    do
	vcftools --vcf $1.orig.$f.recode.vcf --diff $1.bcbio.M.$f.recode.vcf --diff-site --out $1.bcbio.$f
        cat $1.bcbio.$f.diff.sites_in_files | awk '{if ($4==1) print $1"\t"$2}'  > $1.bcbio.$f.lostSNPs.pos
	vcftools --vcf $1.orig.$f.recode.vcf --positions $1.bcbio.$f.lostSNPs.pos --recode --recode-INFO-all --out $1.bcbio.lost.$f
	bedtools intersect -header -a $1.bcbio.lost.$f.recode.vcf -b ~/work/cheo.variants/data/BED/omim.orphanet/omim.bed > $1.bcbio.lost.$f.in_omim.vcf
    done;
}


function new_variants
{
    cat $1.bcbio.SNP.dp10.diff.sites_in_files | awk '{if ($4==2) print $1"\t"$3}' > $1.bcbio.new.dp10.pos
    vcftools --vcf $1.bcbio.M.SNP.dp10.recode.vcf --remove-filtered-all --positions $1.bcbio.new.dp10.pos --recode --recode-INFO-all --out $1.bcbio.new

#jacek db 4644904 variants

    vcftools --vcf $1.bcbio.new.recode.vcf --exclude-positions jacekdb.pos --recode --recode-INFO-all --out $1.bcbio.new.notinjacekdb

    bedtools intersect -header -a $1.bcbio.new.notinjacekdb.recode.vcf -b ~/work/cheo.variants/data/BED/omim.orphanet/omim.bed > $1.bcbio.new.in_omim.vcf
}

#prepare $1
#lost_variants $1
new_variants $1