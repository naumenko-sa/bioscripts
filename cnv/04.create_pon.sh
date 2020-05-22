#!/bin/bash
# $1 = bam.hdf5
# $2 = target.gcannotated.tsv

bname=`basename $1 .hdf5`

unset JAVA_HOME && \
export PATH=/bcbio/anaconda/bin:"$PATH" && \
gatk --java-options '-Xms500m -Xmx131859m -XX:+UseSerialGC -Djava.io.tmpdir=.' \
CreateReadCountPanelOfNormals \
-O $bname.pon.hdf5 \
--annotated-intervals $2 \
-I $1 \
--maximum-zeros-in-sample-percentage 100
