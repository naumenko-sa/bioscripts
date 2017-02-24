#!/bin/bash

# prepares a run of multiples families to run variant calling, one family may have severa samples
# $1 - a file from DCC Re-analysis Queue Final google spreadsheet, columns 2,5,7:
# mapped_sample		family_id	bam
# mapped sample id is unique, but the real family id is project_cohort, because a project might have many family names
# creates one project per family
# run with bcbio.prepare_families.sh table.txt &> file.log to track failed bams

prepare_sample()
{
    local sample=$1

    mkdir -p ${sample}/input
    mkdir ${sample}/work

    cp ~/bioscripts/bcbio.sample_sheet_header.csv $sample.csv

    file=`cat $sample.txt | awk '{print $5}'`
    
    ln -s $file ${sample}/input/${sample}.bam
    
    echo $sample","$sample",,,," >> $sample.csv
    
    bcbio_nextgen.py -w template ~/bioscripts/bcbio.alignment.template.yaml ${sample}.csv ${sample}/input/${sample}.bam
    
    rm $sample.csv
}

cat $1 | awk '{print $2}' >  projects.txt

#create a file with path to bam file for each sample
for sample in `cat projects.txt`
do
    cat $1 | grep $sample > $sample.txt
    prepare_sample $sample
    rm $sample.txt
done
