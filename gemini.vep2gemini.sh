#/bin/bash

#vep2gemini - loads vep annotated vcf file to gemini database 
#based on  bcbio.log

#PBS -l walltime=23:00:00,nodes=1:ppn=16
#PBS -joe .
#PBS -d .
#PBS -l vmem=50g,mem=50g

if [ -z $vcf ]
then
    vcf=$1
fi

bname=`echo $vcf | sed s/.vcf.gz//`

#--skip-cadd if no cadd
#remove GL chromosomes - sometimes they cause problems
mv $vcf $bname.tmp.vcf.gz
gunzip -c $bname.tmp.vcf.gz | grep "^#" > $bname.vcf
gunzip -c $bname.tmp.vcf.gz | grep -v "^#" | grep -v -i "^GL00" >> $bname.vcf
bgzip $bname.vcf
tabix $bname.vcf.gz
rm $bname.tmp.vcf.gz

#remove --passonly to load all variants
#add   --skip-cadd to remove cadd
gemini load --passonly --skip-gerp-bp -v $vcf -t VEP --cores 16 --tempdir . $bname.db



