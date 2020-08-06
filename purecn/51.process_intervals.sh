#!/bin/bash -l

#SBATCH --job-name=purecn
#SBATCH --mem=20G
#SBATCH --export=ALL
#SBATCH -t 7-50:00
#SBATCH -p core -n 1

date

. /projects/ngs/local/users/kmhr378/2020-07-09_BS_benchmark/.bash_profile
bcbio=/projects/ngs/local/users/kmhr378/2020-07-09_BS_benchmark/bcbio

### 5.1 Process intervals file
# $1 = panel.bed
# $2 = GCA_000001405.15_GRCh38_no_alt_analysis_set_100.bw
PURECN=$bcbio/anaconda/envs/r36/lib/R/library/PureCN/extdata

export PATH=$bcbio/anaconda/envs/r36/bin:$PATH

which Rscript

bname=`basename $1 .bed`

Rscript \
$PURECN/IntervalFile.R \
--infile $1 \
--fasta $bcbio/genomes/Hsapiens/hg38/seq/hg38.fa \
--outfile $bname.txt \
--offtarget \
--genome hg38 \
--export panel.optimized.bed \
--mappability /projects/ngs/local/users/kmhr378/2020-06-30_DEV1534_purecn/02_reference/GCA_000001405.15_GRCh38_no_alt_analysis_set_100.bw

date
