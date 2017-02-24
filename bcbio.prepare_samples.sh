#!/bin/bash

# prepares a run of multiples samples to generate bam files
# $1 - a file from DCC Re-analysis Queue Final google spreadsheet, columns 2,7:
# mapped_sample		bam
# mapped sample id is unique
# creates one project per sample
# for the preparation of variant calling on the per family basis use another script and template
# run with bcbio.prepare_samples.sh table.txt &> file.log to track failed bams

# debug - after refactoring!

prepare_sample()
{
    local sample=$1
    local bam=$2

    mkdir -p ${sample}/input
    mkdir ${sample}/work

    cp ~/bioscripts/bcbio.sample_sheet_header.csv $sample.csv

    ln -s $bam ${sample}/input/${sample}.bam
    
    #no batch for single sample
    echo $sample","$sample",,,," >> $sample.csv
    
    bcbio_nextgen.py -w template ~/bioscripts/bcbio.alignment.template.yaml ${sample}.csv ${sample}/input/${sample}.bam
    
    rm $sample.csv
}


while read sample bam
do
    prepare_sample $sample $bam
done < $1
