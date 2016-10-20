#!/bin/bash

#validate vcf file against genome in a bottle calls for NA12878

tabix $1

cat /hpf/largeprojects/ccmbio/naumenko/NA12878-exome-eval/input/GiaB_v2_19_regions.bed | grep -v ^track | grep -v ^browser | \
   grep -v "^#" | /home/naumenko/work/tools/bcbio/anaconda/bin/py -x 'bcbio.variation.bedutils.remove_bad(x)' | sort -V -k1,1 -k2,2n > GiaB_v2_19_regions.bed

cat GiaB_v2_19_regions.bed | /home/naumenko/work/tools/bcbio/anaconda/bin/pbgzip -n 5  -c > GiaB_v2_19_regions.bed.gz

/home/naumenko/work/tools/bcbio/anaconda/bin/tabix -f -p bed GiaB_v2_19_regions.bed.gz

bedtools intersect -nonamecheck -a NA12878-sort-callable_sample.bed -b GiaB_v2_19_regions.bed > NA12878-sort-callable_sample-NA12878-wrm.bed

#uses PASS variants only
export PATH=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/anaconda/bin:$PATH &&  \
   export RTG_JAVA_OPTS='-Xms750m' && export RTG_MEM=9100m && \
   rtg vcfeval --threads 5 -b /hpf/largeprojects/ccmbio/naumenko/NA12878-exome-eval/input/GiaB_v2_19.vcf.gz \
   --bed-regions NA12878-sort-callable_sample-NA12878-wrm.bed \
   -c $1 \
   -t /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/rtg/GRCh37.sdf \
   -o rtg --vcf-score-field='GQ' 
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

