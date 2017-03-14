#!/bin/bash

if [ -z $family];
then
    family=$1
fi

cd $family
mkdir input;

mv *-*/*.gz input/

cp ~/bioscripts/bcbio.sample_sheet_header.csv $family.csv

#rename files
cd input
for f in *.fastq.gz;
do 
    mv $f `echo $f | sed s/"_001.fastq.gz"/.fq.gz/ | sed s/"S._R"//`;
done;
cd ..

ls -1 input/*1.fq.gz | awk -F '/' '{print $2}' | awk -F '_' -v project=$family '{print $1","$1","project",,,"}' >> $family.csv

cd ..

bcbio_nextgen.py -w template ~/bioscripts/bcbio.templates.exome.yaml ${family}/${family}.csv ${family}/input/*

#cd $family/work

#qsub ~/bioscripts/bcbio.pbs -v project=$family

#cd ..
