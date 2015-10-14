#/bin/bash
date
export PATH=/mnt/lustre/tools/rsem-1.2.19:$PATH
#rsem-prepare-reference --bowtie gam8.5.fasta gam8.5
rsem-calculate-expression -p 20  --paired-end gam8.5.filt.r1.trim.fq gam8.5.filt.r2.trim.fq gam8.5 gam8.5
date