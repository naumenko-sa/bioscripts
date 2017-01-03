#!/bin/bash

# cleans up after bcbio - when running large cohort only final folder is kept, and only ensemble gemini database: 2-3G per family
# prepares tables for report generation
# creates report

# parameters:
# family = [family_id] (=folder_name)

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

function cleanup
{

    mv ${family}/final/2016*/* $family
    rm -rf ${family}/final/
    rm -rf ${family}/work/
    rm -rf ${family}/input/

    rm ${family}/${family}-freebayes.db
    rm ${family}/${family}-gatk-haplotype.db
    rm ${family}/${family}-samtools.db
    rm ${family}/${family}-platypus.db
}

function prepare_for_report
{
    cd $family

    gemini.gemini2txt.sh ${family}-ensemble.db 

    gemini.decompose.sh ${family}-freebayes.vcf.gz
    vcf.freebayes.getAO.sh ${family}-freebayes.decomposed.vcf.gz

    gemini.decompose.sh ${family}-gatk-haplotype.vcf.gz
    vcf.gatk.get_depth.sh ${family}-gatk-haplotype.decomposed.vcf.gz

    gemini.decompose.sh ${family}-platypus.vcf.gz
    vcf.platypus.getNV.sh ${family}-platypus.decomposed.vcf.gz

    cd ..
}

cleanup
prepare_for_report
