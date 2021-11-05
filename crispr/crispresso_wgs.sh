#!/bin/bash -l

#SBATCH --job-name=crispresso
#SBATCH --mem=20G
#SBATCH --export=ALL
#SBATCH -t 7-50:00
#SBATCH -p core -n 1

date
# $1 = sample
ml bcbio-nextgen/1.2.8

# finally, baits should be padded 100 bp first and then sliced into chunks
# because crispresso is for amplicons (start at the same position),
# when using panels, the window should be small, since it is cuttings reads on the window
# if the window is large, there is fewer reads covering it - no coverage
# It could be a situation like  [---- bait1 --- ] [ editing site] [ --- bait2 ---------]
# so if supplying only baits, the editing site is missing, so expand, merge, and slice into
#                               [ --- slice1 ] [ slice 2] [slide 3] [ slide4 ][ slice5 ]

# bedtools makewindows -b baits.mm10.sorted.padded100bp.merged.bed -w 40 > baits.mm10.sorted.padded100bp.merged.split_window40.bed
# cat baits.mm10.sorted.padded100bp.merged.split_window40.bed | awk '{print $0"\tbait"NR"\tNA\tNA\tNA"}' > baits.mm10.sorted.padded100bp.merged.split_window40.tsv

# note the visualization parameter - you could miss offtargets with the default 0.2%

singularity run \
-e crispresso2_latest.sif CRISPRessoWGS \
-b ${1}-ready.bam \
-f baits.mm10.sorted.padded100bp.merged.split_window40.tsv \
-r /projects/ngs/reference/UpdateGenomesBcbio/Mmusculus/mm10_gene/seq/mm10_gene.fa \
-n ${1}_wgs_mode_window_center_6 \
--quantification_window_center -6 \
--debug \
--min_frequency_alleles_around_cut_to_plot 0.0001


date
