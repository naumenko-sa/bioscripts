###############################################################################
# Filter TCAG CNV report
###############################################################################
library(tidyverse)

args = commandArgs(trailingOnly = T)
input_table <- read_tsv(args[1], col_types = cols(.default = "c"))
panel <- read_csv(args[2])

# gene_symbol = "TPTE|BAGE4"
# report CNV is at least one gene overlaps
filtered_table <- NULL
for(i in 1:nrow(input_table)){
    genes <- input_table[i, "gene_symbol"]
    if (!is.na(genes)){
	genes_array <- str_split(genes, "\\|")[[1]]
	shared_genes <- intersect(genes_array, panel$external_gene_name)
	if (length(shared_genes) > 0){
	    filtered_table <- bind_rows(filtered_table, slice(input_table, i))
	}
    }
}

write_excel_csv(filtered_table, args[3])
