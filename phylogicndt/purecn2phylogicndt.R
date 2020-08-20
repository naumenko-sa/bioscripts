library(tidyverse)

# local_cn_a1: Minor allele local copy number = M, coverage_loess_genes.csv
# local_cn_a2: Major allele local copy number = C, coverage_loess_genes.csv

# usage: Rscript coverage_loess_genes.csv coverage_loess_variants.csv coverage_loess_loh.csv sample.maf.tsv sample.cn.tsv
# output: sample.maf.tsv sample.cn.tsv for phylgicndt

args = commandArgs(trailingOnly=TRUE)

genes <- read_csv(args[1]) %>% 
    dplyr::select(gene.symbol, M, C)

variants <- read_csv(args[2]) %>% 
    dplyr::filter(FLAGGED == FALSE) %>% drop_na(gene.symbol) %>% 
    mutate(t_ref_count = round(depth* (1-AR)), t_alt_count = round(depth * AR)) %>% 
    dplyr::select(gene.symbol, chr, start, REF, ALT, t_ref_count, t_alt_count) %>% 
    left_join(genes, by = c("gene.symbol" = "gene.symbol")) %>% 
    dplyr::rename(Hugo_Symbol = gene.symbol,
                  Chromosome = chr,
                  Start_position = start,
                  Reference_Allele = REF,
                  Tumor_Seq_Allele2 = ALT,
                  local_cn_a1 = M,
                  local_cn_a2 = C) %>% 
    write_tsv(args[4])

segments <- read_csv(args[3]) %>%
    dplyr::select(chr, start, end, M, C) %>%
    drop_na(M, C) %>%
    dplyr::rename(Chromosome = chr, Start_position = start, End_Position = end, A1_CN = M, A2_CN = C) %>%
    dplyr::mutate(ID = 1:n(), A1_CN = round(A1_CN), A2_CN = round(A2_CN)) %>%
    dplyr::select(ID, everything()) %>%
    write_tsv(args[5])
