#!/bin/bash

#compares 2 vcf files

module load bcftools/1.3
bcftools isec -p _dir $1 $2

echo "Variants:"
echo "Unique to" $1 `cat _dir/0000.vcf | grep -c -v "^#"`
echo "Unique to" $2 `cat _dir/0001.vcf | grep -c -v "^#"`
echo "Shared" `cat _dir/0002.vcf | grep -c -v "^#"`

