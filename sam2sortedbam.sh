#!/bin/bash

bamname=`echo $1 | sed s/sam/bam/`;

samtools view -bS $1 > $bamname
sortedname=`echo $1 | sed s/sam/sorted/`
samtools sort $bamname $sortedname
samtools index $sortedname.bam