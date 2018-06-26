#!/bin/bash

# compares 2 vcf files, needs bcftools in PATH

for arg do
    echo $arg
    bcftools stats $arg | egrep "number of SNPs:" | awk '{print $NF"\tSNPs"}'
    bcftools stats $arg | egrep "number of indels:" | awk '{print $NF"\tindels"}'
done

mkdir _dir
bcftools isec -p _dir $1 $2

echo "Variants:"
echo "Unique to" $1 `cat _dir/0000.vcf | grep -c -v "^#"` "AVR QUAL=" `cat _dir/0000.vcf | grep -v "^#" | awk '{sum+=$6}END{print sum/NR}'`
echo "Unique to" $2 `cat _dir/0001.vcf | grep -c -v "^#"` "AVR QUAL=" `cat _dir/0001.vcf | grep -v "^#" | awk '{sum+=$6}END{print sum/NR}'`
echo "Shared1" `cat _dir/0002.vcf | grep -c -v "^#"` "AVR QUAL=" `cat _dir/0002.vcf | grep -v "^#" | awk '{sum+=$6}END{print sum/NR}'`
echo "Shared2" `cat _dir/0003.vcf | grep -c -v "^#"` "AVR QUAL=" `cat _dir/0003.vcf | grep -v "^#" | awk '{sum+=$6}END{print sum/NR}'`
