#!/bin/bash

bname1=`basename $1 .vcf.gz`
bname2=`basename $2 .vcf.gz`

bcftools isec $1 $2 -p ${bname1}_${bname2} -n =2 -w 1
