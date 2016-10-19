#!/bin/bash

#runs VEP annotation for RNA-seq mutations 
#based on bcbio.log

unset PERL5LIB && export PATH=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/anaconda/bin:$PATH && /home/naumenko/work/tools/bcbio/anaconda/bin/variant_effect_predictor.pl --vcf -o stdout \
    -i $1 --fork 16 --species homo_sapiens --no_stats --cache --offline --dir /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/vep \
    --symbol --numbers --biotype --total_length --canonical --ccds --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,CANONICAL,CCDS,LoF,LoF_filter,LoF_flags \
    --sift b --polyphen b --plugin LoF,human_ancestor_fa:/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/variation/human_ancestor.fa.gz | sed '/^#/! s/;;/;/g' | bgzip -c > \
    `echo $1 | sed s/vcf/vepeffects.vcf/`
