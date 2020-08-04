#!/bin/bash

#SBATCH --job-name=bcbio
#SBATCH --mem=10G
#SBATCH --export=ALL
#SBATCH -t 7-50:00
#SBATCH -p core -n 1

# gathers all hdf5 files in a PON
# $1 = panel.gcannotated.tsv

hdf5_files=""
for f in *.hdf5
do
    hdf5_files="$hdf5_files -I $f"
done

unset JAVA_HOME && \
export PATH=/bcbio/anaconda/bin:"$PATH" && \
gatk --java-options '-Xms500m -Xmx131859m -XX:+UseSerialGC -Djava.io.tmpdir=.' \
CreateReadCountPanelOfNormals \
-O cnv.pon.hdf5 \
--annotated-intervals $1 \
$hdf5_files \
--maximum-zeros-in-sample-percentage 100
