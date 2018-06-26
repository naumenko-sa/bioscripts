#!/bin/bash

# compares 2 vcf files, needs bcftools in PATH

function get_snps
{
    echo $1
    bcftools stats $1 | egrep "number of SNPs:" | awk '{print $NF"\tSNPs"}'
}

function get_indels
{
    echo $1
    bcftools stats $arg | egrep "number of indels:" | awk '{print $NF"\tindels"}'
}

for arg do
    get_snps $arg
    get_indels $arg
done

mkdir _dir
bcftools isec -p _dir $1 $2

cd _dir

for f in *.vcf;do bgzip $f;tabix $f.vcf.gz;done;

cd ..

echo "Unique to" $1 
get_snps _dir/0000.vcf.gz
get_indels _dir/0000.vcf.gz

echo "Unique to" $2
get_snps _dir/0001.vcf.gz
get_indels _dir/0001.vcf.gz

echo "Shared"
get_snps _dir/0002.vcf.gz
get_indels _dir/0002.vcf.gz
