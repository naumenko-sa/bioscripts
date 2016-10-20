#/bin/bash

#vep2gemini - for rna-seq based on  bcbio.log

#PBS -l walltime=2:00:00,nodes=1:ppn=16
#PBS -joe .
#PBS -d .
#PBS -l vmem=50g,mem=50g

set $vcf

bname=`echo $1 | sed s/.vcf.gz//`

#gunzip -c $1 | sed 's/ID=AD,Number=./ID=AD,Number=R/' | vt decompose -s - | awk '{ gsub("./-65", "./."); print $0 }'| bgzip -c > $bname.decompose.vcf.gz

#/home/naumenko/work/tools/bcbio/anaconda/bin/tabix -f -p vcf $bname.decompose.vcf.gz

#unset PERL5LIB && export PATH=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/anaconda/bin:$PATH && /home/naumenko/work/tools/bcbio/anaconda/bin/variant_effect_predictor.pl --vcf -o stdout \
#    -i $bname.decompose.vcf.gz --fork 16 --species homo_sapiens --no_stats --cache --offline --dir /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/vep \
#    --symbol --numbers --biotype --total_length --canonical --ccds --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,CANONICAL,CCDS,LoF,LoF_filter,LoF_flags \
#    --sift b --polyphen b --plugin LoF,human_ancestor_fa:/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/variation/human_ancestor.fa.gz | sed '/^#/! s/;;/;/g' | \
#    bgzip -c > $bname.decompose.vepeffects.vcf.gz

#/home/naumenko/work/tools/bcbio/anaconda/bin/tabix -f -p vcf $bname.decompose.vepeffects.vcf.gz

/home/naumenko/work/tools/bcbio/anaconda/bin/gemini load  --passonly --skip-cadd --skip-gerp-bp  -v $vcf -t VEP --cores 16 --tempdir . $bname.db



