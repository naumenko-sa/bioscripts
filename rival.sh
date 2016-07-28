#!/bin/bash
#PBS -l walltime=50:00:00,nodes=1:ppn=10
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g

####################################################################
###  GATK variant calling best practices using RNA-seq data
###  https://www.broadinstitute.org/gatk/guide/article?id=3891
###  for human genome

###  Developed for RIVAL project at SickKids 
###  Requires bioscripts
####################################################################

echo "START: " `date`

export PICARD=/hpf/tools/centos6/picard-tools/2.0.1/picard.jar
export JAVA=/hpf/tools/centos6/java/1.8.0_65/bin/java
export GATK=/hpf/tools/centos6/gatk/3.5.0/GenomeAnalysisTK.jar
export GENOME=/home/naumenko/work/reference_2pass
export THREADS=20
export REFERENCE=/hpf/largeprojects/ccmbio/arun/Tools/Genomes/hg38/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa
export DBSNP=/hpf/largeprojects/ccmbio/arun/Tools/Genomes/hg38/hg38bundle/dbsnp_144.hg38.vcf.gz
#export REFERENCE=/home/naumenko/work/rocket/reference/hg38.chr20.fasta
export srr=SRR307898
export SNPEFF=/hpf/largeprojects/ccmbio/arun/Tools/SNPEff/snpEff4.2/snpEff.jar
export RAM=10g

function picard_validate_sam
{
    $JAVA -jar $PICARD ValidateSamFile I=$1 MODE=SUMMARY
}


function step1_create_index_n_dictionary
{
    $JAVA -jar $PICARD CreateSequenceDictionary R=$REFERENCE O=`echo $REFERENCE | sed s/fasta/dict/`
    samtools faidx $REFERENCE
}

#first run prefetch SRRNNNN on data1 node to download sra file to ~/ncbi
function step1_get_fastq
{
    fastq-dump -I --split-files $srr
    rm ~/ncbi/public/sra/$srr.sra
}

#Q>=13, len>=50
function step2_trim
{
    >adapters.fasta;
    fastq.trim.sh ${srr}_1.fastq ${srr}_2.fastq adapters.fasta 13 50 $THREADS
    mv ${srr}_1.fastq.trim ${srr}_1.fq
    mv ${srr}_2.fastq.trim ${srr}_2.fq
    rm ${srr}_1.fastq ${srr}_2.fastq
}

#function step3_mapping
#{
    #use separate STAR script
#}

#$! - input.bam, $2 - output.bam
function step4_remove_duplicates
{
    #sam2bam and mark duplicates by picard
    $JAVA -Xmx4g -jar $PICARD AddOrReplaceReadGroups INPUT=$1 OUTPUT=$1.tmp \
	SO=coordinate RGID=id RGLB=library RGPL=ILLUMINA RGPU=machine RGSM=sample

    $JAVA -Xmx4g -jar $PICARD MarkDuplicates INPUT=$1.tmp OUTPUT=$2 \
     CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT METRICS_FILE=output.metrics \
     # it is debated should we really remove them or just mark
     REMOVE_DUPLICATES=true
     
     rm $1.tmp
}
#$1 - input.bam, $2 - output.bam
function step5_split_n_trim
{
    $JAVA -Xmx4g -jar $GATK -T SplitNCigarReads -R $REFERENCE \
   -I $1 -o $2 -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS \
   -et NO_ET -K /home/naumenko/tools/evolgenomicslab_gmail.com.key \
   -fixNDN
}

function step6_base_recalibration
{
    $JAVA -Xmx4g -jar $GATK -T BaseRecalibrator -R $REFERENCE \
	-I $1 -o recalibration_report.grp -et NO_ET -K /home/naumenko/tools/evolgenomicslab_gmail.com.key \
        -knownSites $DBSNP
        
    $JAVA -Xmx4g -jar $GATK -T PrintReads -R $REFERENCE -I $1 \
	-BQSR recalibration_report.grp -o $2 \
	-et NO_ET -K /home/naumenko/tools/evolgenomicslab_gmail.com.key
}

function step7_call_haplotypes
{
    $JAVA -Xmx4g -jar $GATK -T HaplotypeCaller \
       -R $REFERENCE -I $1 -dontUseSoftClippedBases \
       -stand_call_conf 20.0 -stand_emit_conf 20.0 \
        -o $2 \
       -et NO_ET -K /home/naumenko/tools/evolgenomicslab_gmail.com.key

}

function step8_variant_filtration
{
    $JAVA -Xmx4g -jar $GATK -T VariantFiltration \
    -R $REFERENCE -V $1 -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" \
    -filterName QD -filter "QD < 2.0" -o $1.tmp \
    -et NO_ET -K /home/naumenko/tools/evolgenomicslab_gmail.com.key
    
    #filter out variants without PASS in the filter field
    vcftools --vcf $1.tmp --remove-filtered-all --recode --recode-INFO-all --out $2
    rm $1.tmp
}

function step9_annotateWithSnpeff
{
    $JAVA -Xmx4g -jar $SNPEFF GRCh38.82 $1 > $2
}

function more_steps
{
#remove non-coding and non-splicing variants
cat $srr.ann.vcf | egrep "^#|frameshift_variant|missense_variant|splice|start_lost|stop_gained|stop_lost|synonymous_variant" > $srr.coding.vcf;

#filter dbsnp database
cat $srr.coding.vcf | grep -v '^#'  | awk '{print $1"\t"$2}' | sed s/chr// > $srr.coding.list;
while read chr pos;
do
    tabix $GENOME/mgp.v5.merged.snps_all.dbSNP142.vcf.gz $chr:$pos-$pos >> $srr.dbsnp.vcf;
done < $srr.coding.list;

cat $srr.coding.vcf | grep -v '^#' | awk '{print $1"\t"$2}' | sed s/chr// | sort > $srr.coding.sorted;
cat $srr.dbsnp.vcf | grep -v '^#' | awk '{print $1"\t"$2}' | sort > $srr.dbsnp.sorted;

#cat $srr.coding.vcf | grep "^#" > $srr.mutations.vcf;
comm -23 $srr.coding.sorted $srr.dbsnp.sorted | awk '{print "chr"$1"\t"$2}' > $srr.mutations.list;
vcftools --positions $srr.mutations.list --recode --recode-INFO-all --vcf $srr.ann.vcf --out $srr.mutations;
mv $srr.mutations.recode.vcf $srr.mutations.vcf;
}

#step1_create_index_n_dictionary
#step4_remove_duplicates $srr.chr20.bam $srr.chr20.dedup.bam
#step5_split_n_trim $srr.chr20.dedup.bam $srr.chr20.split.bam
#picard_validate_sam $srr.chr20.bam
#step6_base_recalibration $srr.chr20.split.bam $srr.chr20.recalibrated.bam
#step7_call_haplotypes $srr.chr20.recalibrated.bam $srr.chr20.vcf
#step8_variant_filtration $srr.chr20.vcf $srr.chr20
step9_annotateWithSnpeff $srr.chr20.recode.vcf $srr.chr20.annotated.vcf
#annotateWithSnpeff NA12878.chr20.vcf NA12878.chr20.annotated.vcf
echo "END: " `date`