#!/bin/bash

cat /hpf/largeprojects/ccmbio/naumenko/NA12878-exome-eval/input/GiaB_v2_19_regions.bed | grep -v ^track | grep -v ^browser | \
   grep -v ^# \
   | /home/naumenko/work/tools/bcbio/anaconda/bin/py -x 'bcbio.variation.bedutils.remove_bad(x)' | sort -V -k1,1 -k2,2n > \
   /hpf/largeprojects/ccmbio/naumenko/NA12878-exome-eval/eval_genap/GiaB_v2_19_regions.bed

cat /hpf/largeprojects/ccmbio/naumenko/NA12878-exome-eval/eval_genap/GiaB_v2_19_regions.bed | /home/naumenko/work/tools/bcbio/anaconda/bin/pbgzip -n 5  -c > \
   /hpf/largeprojects/ccmbio/naumenko/NA12878-exome-eval/eval_genap/GiaB_v2_19_regions.bed.gz

/home/naumenko/work/tools/bcbio/anaconda/bin/tabix -f -p bed /hpf/largeprojects/ccmbio/naumenko/NA12878-exome-eval/eval_genap/GiaB_v2_19_regions.bed.gz

bedtools intersect -nonamecheck -a NA12878-sort-callable_sample.bed -b GiaB_v2_19_regions.bed > \
   /hpf/largeprojects/ccmbio/naumenko/NA12878-exome-eval/eval_genap/NA12878-sort-callable_sample-NA12878-wrm.bed

#uses PASS variants only
export PATH=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/anaconda/bin:$PATH &&  \
   export RTG_JAVA_OPTS='-Xms750m' && export RTG_MEM=9100m && \
   rtg vcfeval --threads 5 -b /hpf/largeprojects/ccmbio/naumenko/NA12878-exome-eval/input/GiaB_v2_19.vcf.gz \
   --bed-regions /hpf/largeprojects/ccmbio/naumenko/NA12878-exome-eval/eval_genap/NA12878-sort-callable_sample-NA12878-wrm.bed \
   -c $1 \
   -t /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/rtg/GRCh37.sdf \
   -o /hpf/largeprojects/ccmbio/naumenko/NA12878-exome-eval/eval_genap/rtg --vcf-score-field='GQ'

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

