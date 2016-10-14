#ask for a complete database from your bioinformatics provider - no medically-relevant filtration

library(dplyr)
setwd("~/Desktop")
file="MW_159_z1657-gatk-haplotype.txt"

full_set <- read.delim(file, stringsAsFactors=FALSE)

restricted = filter(full_set,is_coding==1 & is_exonic==1 & impact!='synonymous_variant' & aaf_exac_all<0.01)

write.table(restricted,paste0(file,".filtered"),quote=F,row.names=F)


