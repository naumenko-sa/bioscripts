#!/bin/bash

# prepares a run of multiples samples to generate bam files
# $1 - file from DCC Re-analysis Queue Final google spreadsheet, columns 2-4
# creates one project per sample
# for the preparation of variant calling on the per family basis use another script and template
# run with cheo.prepare_projects.sh table.txt &> file.log for the case of failed projects

prepare_sample()
{
    local sample=$1

    mkdir -p ${sample}/input
    mkdir ${sample}/work

    cp ~/bioscripts/cheo.prepare.header.txt $sample.csv

    file=`cat $sample.txt | awk '{print $5}'`
    
    ln -s $file ${sample}/input/${sample}.bam
    
    echo $sample","$sample",,,," >> $sample.csv
    
    bcbio_nextgen.py -w template ~/bioscripts/bcbio.alignment.template.yaml ${sample}.csv ${sample}/input/${sample}.bam
    
    rm $sample.csv
}

cat $1 | awk '{print $1}' >  projects.txt

#create a file with path to bam file for each sample
for sample in `cat projects.txt`
do
    cat $1 | grep $sample > $sample.txt
    prepare_sample $sample
    rm $sample.txt
done
