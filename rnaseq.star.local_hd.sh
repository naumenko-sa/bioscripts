#!/bin/bash
#PBS -l walltime=10:00:00,nodes=1:ppn=20
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g

####################################################################
###  GATK variant calling best practices using RNA-seq data
###  https://www.broadinstitute.org/gatk/guide/article?id=3891
###  for human genome
###  STAR mapping part - we need star mapping for all tools
####################################################################

export GENOME=/home/naumenko/work/rocket/reference/;
export THREADS=20
export srr=SRR307897
#STAR is in the PATH

function step1_genome_generate
{
cd $GENOME
STAR --runMode genomeGenerate --genomeDir $GENOME --genomeFastaFiles hg38.chr20.fasta --runThreadN $THREADS
cd -
}

function local_hd_prepare
{
    cp ${srr}_1.fq ${srr}_2.fq $TMP_DIR
    cd $TMP_DIR
}

function local_hd_cleanup
{
    rm $TMP_DIR/*.fq
    mkdir $CURR_DIR/results
    cp -r $TMP_DIR/* $CURR_DIR/results
    rm $TMP_DIR/*
    cd $CURR_DIR
}

function step2_4_map
{
    mkdir 1pass
    cd 1pass

    STAR --genomeDir $GENOME --readFilesIn ../${srr}_1.fq ../${srr}_2.fq --runThreadN $THREADS

    cd ..

    mkdir 2pass

    #overhang should be read length-1
    cp $GENOME/hg38.chr20.fasta 2pass
    cp 1pass/SJ.out.tab 2pass
    cd 2pass

    STAR --runMode genomeGenerate --genomeDir . \
      --genomeFastaFiles hg38.chr20.fasta \
        --sjdbGTFfile /hpf/largeprojects/ccmbio/arun/Tools/Genomes/hg38/gtf/Homo_sapiens.GRCh38.84.gtf \
	 --runThreadN $THREADS \
          --sjdbOverhang 75 \
	    --sjdbFileChrStartEnd SJ.out.tab;
    cd ..

    rm Aligned.out.sam;

    cd 2pass
    STAR --genomeDir 2pass --readFilesIn ../${srr}_1.fq ../${srr}_2.fq --runThreadN $THREADS
}

echo "START: " `date`

export TMP_DIR=/localhd/`echo $PBS_JOBID | cut -d. -f1 `
export CURR_DIR=`pwd`

local_hd_prepare

step2_4_map

local_hd_cleanup

echo "END: " `date`