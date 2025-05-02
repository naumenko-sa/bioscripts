#!/bin/bash

sample=$1

if [ ! -f $sample.vcf.gz ];then
    
    mkdir -p $sample

    cp chr*/${sample}* $sample

    cd $sample

    # females don't have chrY variants and bcftools concat fails
    ny=`gunzip -c $sample.chrY.vcf.gz | grep -vc "^#"`
    if [ $ny == 0 ];then rm $sample.chrY.vcf.gz $sample.chrY.vcf.gz.tbi;fi;

    bcftools concat `ls  -v *vcf.gz` -Ov | bgzip -c > $sample.vcf.gz
    tabix $sample.vcf.gz

    rm $sample.chr*.vcf.gz
    rm $sample.chr*.vcf.gz.tbi

    cd ..

    mv ${sample}/$sample.* .

    rmdir $sample
fi
