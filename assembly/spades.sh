#!/bin/bash
#SBATCH -t 7-00:00:00
#SBATCH -p shared
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --job-name=spades
#SBATCH --mail-type=NONE


date

ml python/2.7.14-fasrc02
ml centos6/0.0.1-fasrc01
ml ncf/1.0.0-fasrc01
ml SPAdes

#spades.py --pe1-1 short.r1.fq  --pe1-2 short.r2.fq  --pe1-s short.single.fq --mp1-1 mp56.r1.fq --mp1-2 mp56.r2.fq --mp2-1 mp67.r1.fq --mp2-2 mp67.r2.fq --mp3-1 mp810.r1.fq --mp3-2 mp810.r2.fq -o spades_komar -t 36
spades.py \
-1 $1 \
-2 $2 \
-o assembly_clean  \
-t 10

date
