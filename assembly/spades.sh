#!/bin/bash

#PBS -d .
#PBS -l walltime=1000:00:00,nodes=1:ppn=23

date
#spades.py --pe1-1 short.r1.fq  --pe1-2 short.r2.fq  --pe1-s short.single.fq --mp1-1 mp56.r1.fq --mp1-2 mp56.r2.fq --mp2-1 mp67.r1.fq --mp2-2 mp67.r2.fq --mp3-1 mp810.r1.fq --mp3-2 mp810.r2.fq -o spades_komar -t 36
spades.py --pe1-1 $1 --pe1-2 $2 -o 3-1-M.PLIN4.spades  -t 1
date

