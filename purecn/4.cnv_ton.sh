#!/bin/bash -l

#SBATCH --job-name=purecn
#SBATCH --mem=20G
#SBATCH --export=ALL
#SBATCH -t 7-50:00
#SBATCH -p core -n 1

date

module load bcbio-nextgen/latest-devel

# $1 = sample_T.bam
# $2 = panel.interval_list
# $3 = pon.hdf5
# $4 = panel.gcannotated.tsv

which gatk

bcbio=/projects/ngs/reference/UpdateGenomesBcbio

bname=`basename $1 .bam`

# 4.1 CollectReadCounts
gatk --java-options '-Xms500m -Xmx131859m -XX:+UseSerialGC -Djava.io.tmpdir=.' \
CollectReadCounts \
-I $1 \
-L $2 \
--interval-merging-rule OVERLAPPING_ONLY \
-O $bname.counts.hdf5

# 4.2 Denoise read counts
gatk --java-options "-Xmx12g" \
DenoiseReadCounts \
-I $bname.counts.hdf5 \
--count-panel-of-normals $3 \
--standardized-copy-ratios $bname.standardizedCR.tsv \
--denoised-copy-ratios $bname.denoisedCR.tsv \
--annotated-intervals $4

# 4.3 Model segments
gatk --java-options "-Xmx4g" ModelSegments \
--denoised-copy-ratios $bname.denoisedCR.tsv \
--output . \
--output-prefix $bname

# 4.4 Call copy ratio
# #1 = T.cr.seg
gatk CallCopyRatioSegments \
--input $bname.cr.seg \
--output $bname.called.seg

# 4.5 Plot modelled segments
# $1 = T.denoisedCR.tsv
# $2 hg38.dict
gatk PlotModeledSegments \
--denoised-copy-ratios $bname.denoisedCR.tsv \
--segments $bname.modelFinal.seg \
--sequence-dictionary $bcbio/Hsapiens/hg38/seq/hg38.dict \
--minimum-contig-length 46709983 \
--output segment_plots \
--output-prefix $bname

date
