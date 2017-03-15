#/bin/bash

# gemini.vep.refseq.sh - annotate vcf with vep before loading into gemini database
# based on  bcbio.log
# uses hgvs notation with refseq transcript coordinates and no --pick - all effects for a gene
# first you have to download refseq annotation
# ftp://ftp.ensembl.org/pub/current_variation/VEP/homo_sapiens_refseq_vep_87_GRCh37.tar.gz
# to
# bcbio/genomes/Hsapiens/GRCh37/vep

#PBS -l walltime=2:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

if [ -z $vcf ];
then
    vcf=$1
fi

bname=`echo $vcf | sed s/.vcf.gz//`

unset PERL5LIB && export PATH=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/anaconda/bin:$PATH && /home/naumenko/work/tools/bcbio/anaconda/bin/variant_effect_predictor.pl --vcf -o stdout \
    -i $vcf --species homo_sapiens_refseq --no_stats --cache --offline --dir /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/vep --symbol --numbers --biotype --total_length \
    --canonical --gene_phenotype --ccds --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,CANONICAL,CCDS,HGVSc,HGVSp \
    --plugin LoF,human_ancestor_fa:/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/variation/human_ancestor.fa.gz --sift b --polyphen b --hgvs --shift_hgvs 1 \
    --fasta /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa \
    | sed '/^#/! s/;;/;/g' | bgzip -c > $bname.vepeffects_refseq.vcf.gz

tabix $bname.vepeffects_refseq.vcf.gz