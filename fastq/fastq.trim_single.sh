#!/bin/bash

#SBATCH --partition=priority        # Partition (queue) priority
#SBATCH --time=2-00:00              # Runtime in D-HH:MM format, 10:00:00 for hours
#SBATCH --job-name=qc               # Job name
#SBATCH -c 10			    # cores
#SBATCH --mem=20G                   # Memory needed per CPU or --mem-per-cpu
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=NONE             # Type of email notification (BEGIN, END, FAIL, ALL)

# $1 = input.fq.gz (single end)

#if [ $# -ne "5" ]
#then
#    echo "Runs trimmomatic"
#    echo "Usage: trim_single.sh single.fq adapters.fasta quality[0-40] minlen threads"
#    exit 1
#fi

bname=`basename $1 .fq.gz`

java -jar /n/data1/cores/bcbio/naumenko/tools/Trimmomatic-0.39/trimmomatic-0.39.jar SE \
-threads 10 \
-phred33 \
$1 \
$1.trimmed.fq.gz \
CROP:61

#ILLUMINACLIP:$2:2:40:15 LEADING:$3 TRAILING:$3 SLIDINGWINDOW:4:$3 MINLEN:$4

