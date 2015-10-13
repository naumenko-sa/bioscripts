#!/bin/bash

chr=chr2L
#ind=DGRP-026
#java -Xmx2g -jar /home/tools/picard-tools-1.88/CreateSequenceDictionary.jar R=$chr.fa O=$chr.dict
#samtools faidx chr2R.fa
for ind in `cat freeze2.vcf.individuals`; 
do 
    date >> $chr.progress;
    echo $ind >> $chr.progress;

#vcf-sort
#remove-malformed
#cat freeze2.vcf.sorted | grep -v '#' | awk '{if ($4==$5) print $1" "$2" "$4" "5;}' > freeze2.vcf.malformed &
#remove_malformed.pl

    java -Xmx2g -jar /mnt/home/tools/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R $chr.fa -T SelectVariants -o $chr.$ind.vcf -sn $ind --variant freeze2.vcf.sorted.corrected
    java -Xmx2g -jar /mnt/home/tools/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R $chr.fa -T FastaAlternateReferenceMaker -o $chr.$ind.fa --variant $chr.$ind.vcf
    align2fasta.pl $chr.$ind.fa | awk -v sp=$ind '{ if ($0 ~ />/) print ">"sp;else print $0;}' > $chr.$ind.fasta
    rm $chr.$ind.fa
done;
