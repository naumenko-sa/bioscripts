#!/bin/bash
#PBS -l walltime=150:00:00,nodes=1:ppn=30
#PBS -joe .
#PBS -d .
#PBS -l vmem=50g,mem=50g

. /hpf/largeprojects/ccmbio/naumenko/tools/bcbio_1.1.5/.profile115
which bcbio_nextgen.py

hostname

bcbio_path=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio_1.1.5

#GRCh37
reference="hg38"

# sometimes cannot upgrade STAR on data nodes - memory is low, or cannot take much CPUs
# can take more than a day
export PATH=${bcbio_path}/bin:$PATH && STAR \
--genomeDir ${bcbio_path}/genomes/Hsapiens/${reference}/star \
--genomeFastaFiles ${bcbio_path}/genomes/Hsapiens/${reference}/seq/${reference}.fa \
--runThreadN 30 --limitGenomeGenerateRAM 30000000000 --genomeChrBinNbits 14 --runMode genomeGenerate --genomeSAindexNbases 14
