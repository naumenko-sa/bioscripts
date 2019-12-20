#!/bin/bash

# for some 300x mito genomes from WES you need more memory

ref=/n/shared_db/bcbio/biodata/genomes/Hsapiens/GRCh37/seq/GRCh37.fa

bname=`basename $1 .MT.bam`
echo $bname

gatk --java-options "-Xmx20g" Mutect2 -R $ref -I $1 -O $bname.unfiltered.vcf -L MT

gatk FilterMutectCalls -R $ref -V $bname.unfiltered.vcf -O $bname.vcf
bgzip $bname.vcf
tabix $bname.vcf.gz

bcftools view -f PASS $bname.vcf.gz -Oz > $bname.pass.vcf.gz

#snpEff -Xms750m -Xmx4g eff \
#-dataDir /n/shared_db/bcbio/biodata/genomes/Hsapiens/hg38/snpeff -hgvs -noLog -i vcf -o vcf \
#-csvStats $bname.snpeff.csv \
#-s $bname.effects-stats.html GRCh38.86 $bname.vcf.gz > $bname.effects.vcf

#snpEff -Xmx4g eff \
#-dataDir /n/shared_db/bcbio/biodata/genomes/Hsapiens/hg38/snpeff GRCh38.86 $1 > $2

