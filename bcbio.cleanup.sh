#!/bin/bash

# cleans up after bcbio - when running large cohort only the final folder is kept
# and only one ensemble gemini database: 2-3G per family
# keeps bam files for new samples
# prepares tables for report generation

# parameters:
# family = [family_id] (=folder_name,main result file should be family-ensemble.db)

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=20g,mem=20g

function cleanup
{

    # better to look for project-summary than hardcode the year
    # keep bam files for new samples
    cd $family
    result_dir=`find final -name project-summary.yaml | sed s/"\/project-summary.yaml"//`
    mv $result_dir/* .
    mv final/*/*.bam .
    mv final/*/*.bai .

    rm -rf final/
    rm -rf work/

    #don't remove input files for new projects
    #rm -rf input/

    # we don't need gemini databases for particular calling algorythms
    rm ${family}-freebayes.db
    rm ${family}-gatk-haplotype.db
    rm ${family}-samtools.db
    rm ${family}-platypus.db

    cd ..
}

function make_report
{
    cd $family

    gemini.gemini2txt.sh ${family}-ensemble.db
    gemini.variant_impacts.sh ${family}-ensemble.db
    gemini.refseq.sh $family

    for f in *.vcf.gz;
    do
	tabix $f;
    done

    #report filtered vcf for import in phenotips
    #note that if there is a multiallelic SNP, with one rare allele and one frequent one, both will be reported in the VCF,
    #and just a rare one in the excel report
    cat ${family}-ensemble.db.txt | cut -f 23,24  | sed 1d | sed s/chr// > ${family}-ensemble.db.txt.positions
    bcftools view -R ${family}-ensemble.db.txt.positions -o ${family}.vcf.gz -O z ${family}-ensemble-annotated-decomposed.vcf.gz

    #decompose first for the old version of bcbio!
    #gemini.decompose.sh ${family}-freebayes.vcf.gz
    vcf.freebayes.getAO.sh ${family}-freebayes-annotated-decomposed.vcf.gz

    #gemini.decompose.sh ${family}-gatk-haplotype.vcf.gz
    vcf.gatk.get_depth.sh ${family}-gatk-haplotype-annotated-decomposed.vcf.gz

    #gemini.decompose.sh ${family}-platypus.vcf.gz
    vcf.platypus.getNV.sh ${family}-platypus-annotated-decomposed.vcf.gz

    module load R
    Rscript ~/bioscripts/gemini.gemini2report.R $family

    cd ..
}

if [ -z $family ]
then
    family=$1
fi

echo $family
cleanup
make_report
