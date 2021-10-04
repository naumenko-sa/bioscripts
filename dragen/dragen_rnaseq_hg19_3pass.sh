#!/bin/bash

# build hash

#dragen --build-hash-table true \
#--ht-build-rna-hashtable true \
#--output-directory /staging/reference/genomes/Hsapiens/hg19 \
#--ht-reference /staging/reference/genomes/Hsapiens/hg19/hg19.fa \
#--ht-alt-liftover /opt/edico/liftover/hg19_alt_liftover.sam

date

sample=$1

mkdir -p dragen_${sample}/pass1

project_dir=/staging/2021-10-03_fusions_3pass

# pass 1
if [ ! -f ${project_dir}/dragen_${sample}/pass1/$sample.SJ.out.tab ]; then

dragen \
--output-directory /staging/2021-10-03_fusions_3pass/dragen_${sample}/pass1 \
-r /staging/reference/genomes/Hsapiens/hg19 \
--output-file-prefix $sample \
-1 /staging/2021-10-03_fusions_3pass/input/${sample}_1.fq.gz \
-2 /staging/2021-10-03_fusions_3pass/input/${sample}_2.fq.gz \
--enable-rna true \
--annotation-file /staging/reference/genomes/Hsapiens/hg19/ref-transcripts.gtf \
--enable-rna-quantification false \
--RGID $sample \
--RGSM $sample \
--intermediate-results-dir /staging/tmp \
--enable-rna-gene-fusion false

fi

# pass 2

mkdir -p dragen_${sample}/pass2

if [ ! -f ${project_dir}/dragen_${sample}/pass2/$sample.bam ]; then

dragen \
--output-directory /staging/2021-10-03_fusions_3pass/dragen_${sample}/pass2 \
-r /staging/reference/genomes/Hsapiens/hg19 \
--output-file-prefix $sample \
-1 /staging/2021-10-03_fusions_3pass/input/${sample}_1.fq.gz \
-2 /staging/2021-10-03_fusions_3pass/input/${sample}_2.fq.gz \
--enable-rna true \
--annotation-file /staging/2021-10-03_fusions_3pass/dragen_${sample}/pass1/$sample.SJ.out.tab \
--enable-rna-quantification false \
--RGID $sample \
--RGSM $sample \
--intermediate-results-dir /staging/tmp \
--enable-rna-gene-fusion false

fi

# pass 3
# no read groups

mkdir -p dragen_${sample}/pass3

dragen \
--output-directory /staging/2021-10-03_fusions_3pass/dragen_${sample}/pass3 \
-r /staging/reference/genomes/Hsapiens/hg19 \
--output-file-prefix $sample \
-b ${project_dir}/dragen_${sample}/pass2/$sample.bam \
--enable-rna true \
--annotation-file /staging/reference/genomes/Hsapiens/hg19/ref-transcripts.gtf \
--enable-rna-quantification false \
--rna-quantification-library-type A \
--rna-quantification-gc-bias true \
--intermediate-results-dir /staging/tmp \
--enable-rna-gene-fusion true

date
