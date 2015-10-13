#!/bin/bash

if [ $# -lt 1 ];then
    echo "Extract adapters from fastqc.zip";
    echo "Usage: fastqc_get_adapters.sh file.fastqc.zip";
fi;

bname=`echo $1 | sed s/.zip//`;

unzip $1

cd $bname

cat fastqc_data.txt | awk 'BEGIN{FL=0}{if($0~"Overrepresented") {FL=1};if ((FL==1) && ($0 ~ "Sequence")) FL=2; if (FL==2 && $0 ~ "END_MODULE") FL=3; if (FL==2)print $0}' | grep -v '#' | awk '{print ">"NR"\n"$1}'> ../$bname.adapters.fasta
cd ..

rm -r $bname