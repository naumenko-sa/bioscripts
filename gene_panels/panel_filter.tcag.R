###############################################################################
# Filter TCAG small variant report
###############################################################################
library(tidyverse)
library(writexl)

args = commandArgs(trailingOnly = T)
input_table <- read_tsv(args[1])
panel <- read_csv(args[2])
filtered_table <- input_table %>% filter(FILTER=="PASS", gene_symbol %in% panel$external_gene_name,
                                         !typeseq_priority %in% c("intronic"))
write_xlsx(filtered_table, args[3])
