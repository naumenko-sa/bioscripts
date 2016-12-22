#!/bin/bash

#PBS -l walltime=10:00:00,nodes=1:ppn=10
#PBS -joe .
#PBS -d .
#PBS -l vmem=50g

#mofidied from bcbio rna-seq log and STAR manual

/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/galaxy/../anaconda/bin/STAR --genomeDir /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/star/ \
    --readFilesIn ${project}_1.fq.gz ${project}_2.fq.gz \
    --twopassMode Basic \
     --runThreadN $threads \
     --outFileNamePrefix $project \
     --outReadsUnmapped Fastx \
     --outFilterMultimapNmax 10 \
     --outStd SAM  \
     --outSAMunmapped Within \
     --outSAMattributes NH HI NM MD AS  \
     --sjdbGTFfile /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/rnaseq/ref-transcripts.gtf  \
     --sjdbOverhang $len  --readFilesCommand zcat  --outSAMattrRGline ID:$project PL:illumina PU:$project SM:ID  --outSAMstrandField intronMotif  --quantMode TranscriptomeSAM  \
     | /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/galaxy/../anaconda/bin/samtools sort -@ 5 -m 1G  -T . \
     -o $project.bam /dev/stdin


