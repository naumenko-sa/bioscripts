#!/bin/bash

#SBATCH --partition=priority        # Partition (queue)
#SBATCH --time=5:00:00              # Runtime in D-HH:MM format
#SBATCH --job-name=bcbio           # Job name
#SBATCH -c 1			    # cores
#SBATCH --mem=10G                   # Memory needed per CPU or --mem
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=NONE             # Type of email notification (BEGIN, END, FAIL, ALL)

# gathers all hdf5 files in a PON
# $1 = panel.gcannotated.tsv

. .profile

hdf5_files=""
for f in `cat hdf5.list`
do
    hdf5_files="$hdf5_files -I $f"
done

unset JAVA_HOME && \
gatk --java-options '-Xms500m -Xmx131859m -XX:+UseSerialGC -Djava.io.tmpdir=.' \
CreateReadCountPanelOfNormals \
-O cnv.pon.hdf5 \
--annotated-intervals $1 \
$hdf5_files \
--maximum-zeros-in-sample-percentage 100
