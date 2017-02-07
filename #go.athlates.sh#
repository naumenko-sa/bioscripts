#!/bin/bash -v

#athlates HLA typer

APATH=/mnt/lustre/tools/Athlates_11012013
MPATH=/mnt/lustre/snaumenko1/project_hla/scripts

#eft.fq right.fq hla.clean.fasta

#gobowtie2.sh hla.clean.fasta 1_shuffled.fq 2_shuffled.fq 100 1000 10 fr #result moved to test.sam

samtools view -b -L $MPATH/hla.A.bed $1.sorted.bam > A.bam
samtools view -h -o A.sam A.bam
sam2sortedbam.sh A.sam

samtools view -b -L $MPATH/hla.non-A.bed $1.sorted.bam > nonA.bam
samtools view -h -o nonA.sam nonA.bam
sam2sortedbam.sh nonA.sam

#typing -bam A.sorted.bam -exlbam non-A.sorted.bam -msa $APATH/db/msa/A_nuc.txt -o test_a > test_a.log.txt
typing -bam A.sorted.bam -msa $APATH/db/msa/A_nuc.txt -o test_a > test_a.log.txt