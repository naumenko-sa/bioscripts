# calculate MAPD
# http://tools.thermofisher.com/content/sfs/brochures/mapd_snp6_whitepaper.pdf
# usage Rscript get_mapd.R dnacopy.seg

library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

segments <- read_tsv(args[1])

prev_chr <- "0"
prev_seg_mean <- 0
delta <- c()
for (i in seq(1, nrow(segments))){
    chr <- segments[i, "chrom"] %>% as.character()
    seg_mean <- segments[i, "seg.mean"] %>% as.double()
    
    if (chr == prev_chr){
        d <- abs(seg_mean - prev_seg_mean)
        delta <- c(delta, d)
    }
    
    prev_chr <- chr
    prev_seg_mean <- seg_mean
}
print(round(median(delta), 2))
