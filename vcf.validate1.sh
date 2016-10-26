#!/bin/bash

#validate vcf file against genome in a bottle calls for NA12878 - whole genome version

#tabix $1

#uses PASS variants only
#export PATH=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/anaconda/bin:$PATH &&  \
#   export RTG_JAVA_OPTS='-Xms750m' && export RTG_MEM=9100m && \
#   rtg vcfeval --threads 5 -b $2 \
#   --bed-regions $3 \
#   -c $1 \
#   -t /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/rtg/GRCh37.sdf \
#   -o rtg --vcf-score-field='GQ' 
#   --all-records

module load bcftools
for f in {tp-baseline,fp,fn};
do
    echo snp $f `bcftools view --types snps rtg/$f.vcf.gz | grep -vc "^#"` >> $1.stat
    echo indels $f `bcftools view --exclude-types snps rtg/$f.vcf.gz | grep -vc "^#"` >> $1.stat
done

#vcfanno /hpf/largeprojects/ccmbio/naumenko/NA12878-exome-eval/work/validate/NA12878-2/freebayes/rtg/tx/tmpO5vSAo/tp-baseline-context.toml \
#    /hpf/largeprojects/ccmbio/naumenko/NA12878-exome-eval/work/validate/NA12878-2/freebayes/rtg/tp-baseline.vcf.gz | bgzip -c > \
#    /hpf/largeprojects/ccmbio/naumenko/NA12878-exome-eval/work/validate/NA12878-2/freebayes/rtg/tx/tmpO5vSAo/tp-baseline-context.vcf.gz

#/home/naumenko/work/tools/bcbio/anaconda/bin/tabix -f -p vcf /hpf/largeprojects/ccmbio/naumenko/NA12878-exome-eval/work/validate/NA12878-2/freebayes/rtg/tx/tmpArtZh4/tp-baseline-context.vcf.gz

#vcfanno /hpf/largeprojects/ccmbio/naumenko/NA12878-exome-eval/work/validate/NA12878-2/freebayes/rtg/tx/tmpSSi2pI/fp-context.toml \
#    /hpf/largeprojects/ccmbio/naumenko/NA12878-exome-eval/work/validate/NA12878-2/freebayes/rtg/fp.vcf.gz | bgzip -c > \
#    /hpf/largeprojects/ccmbio/naumenko/NA12878-exome-eval/work/validate/NA12878-2/freebayes/rtg/tx/tmpSSi2pI/fp-context.vcf.gz

#/home/naumenko/work/tools/bcbio/anaconda/bin/tabix -f -p vcf /hpf/largeprojects/ccmbio/naumenko/NA12878-exome-eval/work/validate/NA12878-2/freebayes/rtg/tx/tmpmafMNP/fp-context.vcf.gz

#vcfanno /hpf/largeprojects/ccmbio/naumenko/NA12878-exome-eval/work/validate/NA12878-2/freebayes/rtg/tx/tmpdWS14A/fn-context.toml \
#    /hpf/largeprojects/ccmbio/naumenko/NA12878-exome-eval/work/validate/NA12878-2/freebayes/rtg/fn.vcf.gz | bgzip -c > \
#    /hpf/largeprojects/ccmbio/naumenko/NA12878-exome-eval/work/validate/NA12878-2/freebayes/rtg/tx/tmpdWS14A/fn-context.vcf.gz

#/home/naumenko/work/tools/bcbio/anaconda/bin/tabix -f -p vcf /hpf/largeprojects/ccmbio/naumenko/NA12878-exome-eval/work/validate/NA12878-2/freebayes/rtg/tx/tmpPK7Sze/fn-context.vcf.gz

