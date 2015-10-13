#!/bin/bash
#PBS -d .
#PBS -l walltime=1000:00:00,nodes=node01:ppn=10,mem=520gb

OMP_NUM_THREADS=10
export OMP_NUM_THREADS

date

#cd /mnt/local/snaumenko/gammarus/gam8_1
/mnt/local/snaumenko/gammarus/gam8_1/blastn_mit/assembly
#/mnt/local/snaumenko/gammarus/gam${n}

/home/tools/velvet/velveth kmer$kmer $kmer -fastq -shortPaired gam8_1.filt.fq -long gam8_1.joined.fq 
/home/tools/velvet/velvetg kmer$kmer -read_trkg yes -ins_length 250
oases kmer$kmer -ins_length 250 -unused_reads yes -cov_cutoff 10

#cd kmer$kmer
#cp Log UnusedReads.fa transcripts.fa /home/snaumenko1/project_gamarus/2assembly/gam${n}/
#cd ..
#rm kmer$kmer/*
#rmdir kmer$kmer
#rm gam${n}.fq

#cd /mnt/home/snaumenko1/project_gamarus/3extract_long_isoform

#extract_longest_isoform.sh  kmer$kmer
#../2assembly/gam${n}

date