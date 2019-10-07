###############################################################################
# Filter TCAG small variant report
###############################################################################
library(tidyverse)

args = commandArgs(trailingOnly = T)
input_table <- read_tsv(args[1])
panel <- read_csv(args[2])
filtered_table <- input_table %>% filter(FILTER=="PASS", gene_symbol %in% panel$external_gene_name)
write_excel_csv(filtered_table, args[3])
