#!/bin/bash

#PBS -l walltime=5:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g

#without -csvStats to produce html reports
#csvStats for multiqc

for f in *.vcf;
do
    java -Xmx4G -jar /hpf/tools/centos6/snpEff/4.0/snpEff.jar GRCh37.75 $f -s ${f}.report.html > `echo $f | sed s/vcf/ann.vcf/`;
done;
