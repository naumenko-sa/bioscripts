#!/bin/bash

# prepares a run of multiples families to run variant calling, one family may have severa samples
# $1 - a file from DCC Re-analysis Queue Final google spreadsheet, columns 2,5,7:
# mapped_sample_id	family_id	bam
# mapped sample id is unique, but the real family id is project_cohort, because a project might have many family names
# creates one project per family
# run with bcbio.prepare_families.sh table.txt &> file.log to track failed bams

prepare_family()
{
    local family=$1

    mkdir -p ${family}/input
    mkdir ${family}/work

    cp ~/bioscripts/bcbio.sample_sheet_header.csv $family.csv

    file=`cat $sample.txt | awk '{print $5}'`
    
    while read sample fam bam
    do
	ln -s $bam ${family}/input/${sample}.bam
        echo $sample","$sample","$family",,," >> $family.csv
    done < $family.txt
                
    bcbio_nextgen.py -w template ~/bioscripts/bcbio.templates.exome.yaml $family.csv ${family}/input/*.bam
    
    rm $family.csv
}

cat $1 | awk '{print $2}' | sort | uniq >  families.txt

for family in `cat families.txt`
do
    # not grep because two family names may overlap
    cat $1 | awk -v fam=$family '{if ($2==fam) print $0}' > ${family}.txt
    prepare_family $family
    rm $family.txt
done
