#/bin/bash

vcftools --vcf $1 --remove-filtered-all --remove-indels --recode --out tmp.pass;

cat tmp.pass.recode.vcf | grep -v '^#' | awk '{print $1"\t"$2"\t"$4"\t"$5}' | grep -v "," > $1.mutations;

echo "Total mutations: " `cat $1.mutations | wc -l`;

for context in {CA,CG,CT,TA,TC,TG,GT,GC,GA,AT,AG,AC};
do
    echo $context `cat $1.mutations | awk '{print $3$4}' | grep -c $context`;
done;

rm tmp.pass.recode.vcf;