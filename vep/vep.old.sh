#/bin/bash

# gemini.vcf2vep annotates vcf with vep before loading to gemini database
# based on  bcbio.log
# uses hgvs notation and no --pick = all effects for a gene
# 10h walltime is not enough for genomes or big multisample vcfs

#PBS -l walltime=23:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=30g,mem=30g

if [ -z $vcf ]
then
    vcf=$1
fi

if [ -z $threads ]
then
    threads=5
fi

bname=`basename $vcf .vcf.gz`

#find reference
reference=`readlink -f $(which bcbio_nextgen.py) | sed s/"anaconda\/bin\/bcbio_nextgen.py"/"genomes\/Hsapiens\/hg38"/`
vep_reference=`readlink -f $(which vep) | sed s/"\/vep"//`

echo "Threads:" $threads
echo "Reference:" $reference

#unset PERL5LIB && vep --vcf -o stdout \
#    -i $vcf --fork $threads --species homo_sapiens --no_stats --cache --offline --dir ${reference}/vep --symbol --numbers --biotype --total_length \
#    --canonical --gene_phenotype --ccds --uniprot --domains --regulatory --protein --tsl --appris --af --max_af --af_1kg --af_esp --af_gnomad --pubmed --variant_class \
#    --allele_number \
#    --fasta ${reference}/seq/GRCh37.fa.gz \
#    --plugin LoF,human_ancestor_fa:${reference}/variation/human_ancestor.fa.gz,loftee_path:$vep_reference \
#    --plugin MaxEntScan,/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/anaconda/share/maxentscan-0_2004.04.21-1 \
#    --plugin SpliceRegion --sift b --polyphen b --hgvs --shift_hgvs 1 --merged \
#    | sed '/^#/! s/;;/;/g' | bgzip -c > $bname.vepeffects.vcf.gz

#unset PERL5LIB && vep --vcf -o stdout \
#    -i $vcf --fork $threads --species homo_sapiens --no_stats --cache --offline --dir ${reference}/vep --symbol --numbers --biotype --total_length \
#    --canonical --gene_phenotype --ccds --uniprot --domains --regulatory --protein --tsl --appris --af --max_af --af_1kg --af_esp --af_gnomad --pubmed --variant_class \
#    --allele_number \
#    --fasta ${reference}/seq/GRCh37.fa.gz \
#    --plugin LoF,human_ancestor_fa:${reference}/human_ancestor.fa.gz,loftee_path:$vep_reference \
#    --plugin MaxEntScan,/hpf/largeprojects/ccmbio/naumenko/tools/bcbio_1.1.5/anaconda/share/maxentscan-0_2004.04.21-1 \
#    --plugin SpliceRegion --sift b --polyphen b --hgvsg --hgvs --shift_hgvs 1 --merged \
#    | sed '/^#/! s/;;/;/g' | bgzip -c > $bname.vepeffects.vcf.gz

date

unset PERL5LIB && vep \
-i $vcf \
--vcf \
-o stdout \
--fork $threads \
--species homo_sapiens \
--no_stats \
--cache \
--offline \
--dir ${reference}/vep \
--symbol \
--numbers \
--fasta ${reference}/seq/hg38.fa.gz \
--merged | sed '/^#/! s/;;/;/g' |  bgzip -c > $bname.vepeffects.vcf.gz

date