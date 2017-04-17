#!/bin/bash

# prepares a run of multiples families to run variant calling, one family may have severa samples
# $1 - a file from DCC Re-analysis Queue Final google spreadsheet, columns 2,5,7:
# sample_id	family_id	absolute_path_to_bam_file
# creates one project per family
# 
# run with 
# bcbio.prepare_families.sh table.txt &> file.log to track failed bams
# or
# qsub ~/bioscripts/bcbio.prepare_families.sh -v project_list=table.txt 

# to create table.txt from a directory of bam files names family_sample.yyy.bam
# for f in *.bam;do echo $f | awk -F "." '{print $1"\t"$0}' | awk -F '_' '{print $2"\t"$0}' | awk -v dir=`pwd` '{print $1"\t"$2"\t"dir"/"$4}' >> ~/table.txt;done;

#PBS -l walltime=20:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

prepare_family()
{
    local family=$1

    mkdir -p ${family}/input
    mkdir ${family}/work

    cp ~/bioscripts/bcbio.sample_sheet_header.csv $family.csv

    while read sample fam bam
    do
	ln -s $bam ${family}/input/${sample}.bam
        echo $sample","$sample","$family",,," >> $family.csv
    done < $family.txt
                
    bcbio_nextgen.py -w template ~/bioscripts/bcbio.templates.exome.yaml $family.csv ${family}/input/*.bam
    
    rm $family.csv
}

if [ -z $project_list ];
then
    project_list=$1
fi

cat $project_list | awk '{print $2}' | sort | uniq >  families.txt

for family in `cat families.txt`
do
    # not grep because two family names may overlap
    cat $project_list | awk -v fam=$family '{if ($2==fam) print $0}' > ${family}.txt
    prepare_family $family
    rm $family.txt
done
