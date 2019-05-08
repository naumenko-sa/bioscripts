args <- commandArgs(trailingOnly = T)

#input_report = "182208.wes.2018-11-30.csv"

panel_path <- "~/bioscripts/gene_panels"

input_report <- args[1]

variants <- read.csv(input_report, stringsAsFactors = F)

lupus <- read.csv(paste0(panel_path,"/hiraki_lupus.csv"), stringsAsFactors = F)

lupus.langefeld.gwas <- read.csv(paste0(panel_path,"/hiraki_lupus_langefeld_gwas.csv"), stringsAsFactors = F)

lupus.eastasians.gwas <- read.csv(paste0(panel_path,"/hiraki_lupus_eastasians_gwas.csv"), stringsAsFactors = F)
primary_immunodeficiency <- read.csv(paste0(panel_path,"/genomics_england_primary_immunodeficiency.csv"), stringsAsFactors = F)

periodic_fever <- read.csv(paste0(panel_path, "/genomics_england_periodic_fever_syndromes.csv"), stringsAsFactors = F)

mas <- read.csv(paste0(panel_path,"/hiraki_mas.csv"), stringsAsFactors = F)

recurrent_fever_SK <- read.csv(paste0(panel_path, "/dplm_recurrent_fever_syndrome.csv"), stringsAsFactors = F)

variants$Lupus_panel <- ifelse(variants$Ensembl_gene_id %in% lupus$ensembl_gene_id,"Lupus_panel",NA)
variants$Lupus_langefeld_gwas <- ifelse(variants$Ensembl_gene_id %in% lupus.langefeld.gwas$ensembl_gene_id,"Lupus_langefeld_gwas",NA)
variants$Lupus_eastasians_gwas <- ifelse(variants$Ensembl_gene_id %in% lupus.eastasians.gwas$ensembl_gene_id,"Lupus_eastasians_gwas",NA)
variants$Primary_immunodeficiency_panel <- ifelse(variants$Ensembl_gene_id %in% primary_immunodeficiency$PanelAPP.EnsemblId.GRch37,"Primary_immunodeficiency_panel",NA)
variants$Periodic_fever_panel <- ifelse(variants$Ensembl_gene_id %in% periodic_fever$PanelAPP.EnsemblId.GRch37,"Periodic_fever",NA)
variants$MAS <- ifelse(variants$Ensembl_gene_id %in% mas$ensembl_gene_id,"MAS",NA)
variants$Recurrent_fever_SK <- ifelse(variants$Ensembl_gene_id %in% recurrent_fever_SK$ensembl_gene_id, "Recurrent_fever_SK",NA)

output_report <- sub(".csv", ".immunopanels.csv", input_report)

write.csv(variants, output_report, row.names = F)