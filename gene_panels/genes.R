###############################################################################
# biomart wrappers to get gene, transcript, exons annotations from ENSEMBL
# http://bioconductor.org/packages/release/bioc/html/biomaRt.html
# https://www.bioconductor.org/packages/devel/bioc/vignettes/biomaRt/inst/doc/biomaRt.html
###############################################################################
installation <- function(){
    # lib = "~/R")
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("biomaRt", version = "3.8")
    install.packages("tidyverse")
    #install.packages("bedr")
}

init <- function(){
    library(biomaRt)    
    library(tidyverse)
}

# grch38 by default
# use host=grch37.ensembl.org for grch37 reference
init_mart_human <- function(host = "useast.ensembl.org"){
    mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = host)
    mart <- useDataset(mart, dataset = "hsapiens_gene_ensembl")
    return(mart)
}

tutorial_init_mart_human <- function(){
    library(biomaRt)    
    #library("readr")
    #library(IRanges)
    #library(GenomicRanges)  
    #library(bedr)
    
    listMarts()
    mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org")
    datasets <- listDatasets(mart)
    mart <- useDataset(mart, dataset = "hsapiens_gene_ensembl")
    attributes <- listAttributes(mart)
    filters <- listFilters(mart)
    
    chromosomes <- getBM(attributes = c("chromosome_name"), mart = mart)
    
    genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
                  filters = c("chromosome_name"),
                  values = list("22"),
                  mart = mart)
    
    return(mart)
}

tutorial_explore_attributes_and_filters <- function(){
    mart <- init_mart_human()
    
    attributes <- listAttributes(mart, 
                                 what = c("name", "description", "fullDescription", "page"))
    unique(attributes$page)
    attributes <- listAttributes(mart, "sequences", 
                                 what = c("name", "description", "fullDescription", "page"))
    filters <- listFilters(mart, c("name", "description", "fullDescription"))
}

# for mm10 = grcm38 reference
init_mart_mouse <- function(){
    library(biomaRt)    
    listMarts()
    mart <- useMart(biomart = "ENSEMBL_MART_MOUSE")
    datasets <- listDatasets(mart)
    mart <- useDataset(mart, dataset = "mc57bl6nj_gene_ensembl")
    
    attributes <- listAttributes(mart, what = c("name", "description", "fullDescription", "page"))
    unique(attributes$page)
    filters <- listFilters(mart, c("name", "description", "fullDescription"))
    return(mart)
}

mouse_get_protein_coding_genes <- function(){
    genes_info <- getBM(attributes=c("chromosome_name", "start_position", "end_position", 
                                     "ensembl_exon_id",
                                     "external_gene_name", "mmusculus_homolog_ensembl_gene"),
                        filters = c("biotype", "external_gene_name"),
                        values = list("protein_coding", "DMD"),
                        mart = mart)
    
    #remove transcripts placed on patches
    genes_info <- genes_info[grep('PATCH',genes_info$chromosome_name,invert=T),]
    #remove HSCHR - alleles
    genes_info <- genes_info[grep('HSCHR',genes_info$chromosome_name,invert=T),]
}

get_promoters <- function(){
  mart <- init_mart_human()
  
  attributes <- listAttributes(mart, 
                               what = c("name", "description", "fullDescription", "page"))
  unique(attributes$page)
  attributes <- listAttributes(mart, "sequences", 
                               what = c("name", "description", "fullDescription", "page"))
  filters <- listFilters(mart, c("name", "description", "fullDescription"))
}

tutorial_explore_marts <- function(){
    library(biomaRt)
    listMarts()
 
    mart <- useMart(biomart = "ENSEMBL_MART_SNP")
    datasets <- listDatasets(mart)
    mart <- useDataset(mart, dataset = "hsapiens_snp")
    attributes <- listAttributes(mart, what = c("name", "description", "fullDescription", "page"))
    listFilters(mart, what = c("name", "description", "options"))
    
    mart <- useMart(biomart = "ENSEMBL_MART_FUNCGEN")
    datasets <- listDatasets(mart)
    mart <- useDataset(mart, dataset = "hsapiens_regulatory_feature")
    attributes <- listAttributes(mart, what = c("name", "description", "fullDescription", "page"))
    filters <- listFilters(mart, what = c("name", "description", "options"))
}

gene_descriptions <- function(mart){
    #mart <- init_mart_human()
    ensembl_w_description <- getBM(attributes = c("ensembl_gene_id",
                                                  "external_gene_name",
                                                  "description", "chromosome_name"),
                                   mart = mart)
    write_excel_csv(ensembl_w_description, "ensembl_w_description.csv")
}

# attribute name_1006 is GO_term
# GO_term takes a while for all genes, demo with chrX
tutorial_gene_descriptions_w_go_term <- function(mart){
    ensembl_w_description <- getBM(attributes=c("ensembl_gene_id",
                                               "external_gene_name",
                                               "description",
                                               "name_1006"),
                                  filters = c("chromosome_name"),
                                  values = list("X"),
                                  mart = mart)
    write.csv(ensembl_w_description, file="ensembl_w_description.csv", row.names=F)
}

# test
# swissprot_id = 'P62701'
gene_name_by_uniprotswissprotid <- function(mart, swissprot_id){
    swissprot_id <- "P62701"
    gene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "uniprotswissprot"),
                               filters = c("uniprotswissprot"),
                               values = swissprot_id,
                               mart = mart)
    return(gene$external_gene_name)
}

refseq_transcripts <- function(mart){    
    protein_coding_genes <- getBM(attributes = c("ensembl_gene_id",
                                                 "refseq_mrna",
                                                 "external_gene_name"),
                                  filters = c("biotype"),
                                  values = list("protein_coding"),
                                  mart = mart)
    colnames(protein_coding_genes) <- c("Ensembl_gene_id", "Refseq_transcript_id", "external_gene_name")
    protein_coding_genes <- protein_coding_genes[protein_coding_genes$Refseq_transcript_id!='',]
    write.csv(protein_coding_genes, "refseq_transcripts.csv", row.names = F)
}

# writes a list of external_gene_names to protein_codin_genes.list
get_protein_coding_genes <- function(mart){
    protein_coding_genes <- getBM(attributes = c("ensembl_gene_id",
                                              "ensembl_transcript_id",
                                              "external_gene_name"),
                                 filters = c("biotype"),
                                 values = list("protein_coding"),
                                 mart = mart)
    colnames(protein_coding_genes) = c("Ensembl_gene_id", "Ensembl_transcript_id",
                                       "external_gene_name")                             
    write.csv(protein_coding_genes, "genes.transcripts.csv", row.names = F)
    #colnames(protein_coding_genes)[2] = 'gene_name'
    #gene names might be not unique - polymorphic regions like NCR3 gene, or bugs like CLN3
    #write.table(unique(sort(protein_coding_genes[,1])),
    #            file="protein_coding_genes.list",
    #            quote=F,row.names=F,col.names=F)
    
    #EXAMPLES:
    #- chromosome_name may be an attribute and a filter - use list for values!
    #filterOptions('biotype',mart)
    #attributePages(mart)
    
    protein_coding_genes <- getBM(attributes = c("ensembl_gene_id",
                                                 "external_gene_name",
                                                 "chromosome_name"),
                                 filters = c("biotype", "chromosome_name"),
                                 values = list("protein_coding",22),
                                 mart = mart)
}

# ensemble_gene_id is a mouse strain specific ID
# for protein coding genes
get_gene_descriptions.mouse <- function(mart){
    mart_mouse <- init_mart_mouse()
    ensembl_w_description <- getBM(attributes = c("mmusculus_homolog_ensembl_gene",
                                               "external_gene_name",
                                               "description"),
                                  filters = c("biotype"), 
                                  values = list("protein_coding"),
                                  
                                  mart=mart_mouse)
    colnames(ensembl_w_description)[1]="ensembl_gene_id"
    ensembl_w_description = ensembl_w_description[ensembl_w_description$ensembl_gene_id != "",]
    
    #in the result there are some duplicates, i.e. ENSMUSG00000074254
    #View(ensembl_w_description[duplicated(ensembl_w_description$ensembl_gene_id),])
    ensembl_w_description = ensembl_w_description[!duplicated(ensembl_w_description$ensembl_gene_id),]
    write.csv(ensembl_w_description,file="ensembl_w_description.mouse.csv",row.names=F)
}

# get coding and noncoding transcripts with their ensembl, refseq, ucsc IDs
get_ensembl_refseq_transcript_ids <- function(mart){
    transcripts <- getBM(attributes = c("ensembl_transcript_id", "refseq_mrna", "refseq_ncrna",
                                        "ucsc", "external_gene_name", "external_transcript_name"),
                               mart = mart,
                               filters = "chromosome_name",
                               values = "X")
                               
    ensembl_refseq <- transcripts[transcripts$refseq_mrna!='',]
    write.csv(ensembl_refseq,file="ensembl_refseq.csv", quote=F, row.names=F)
    
    return(ensembl_refseq)
    
    # more info
    # external_transcript_name
}

# coordinates of protein coding genes (all genes), 
# no duplicate record!
protein_coding_genes_bed <- function(mart){
    genes_info <- getBM(attributes=c("chromosome_name", "start_position", "end_position",
                                     "external_gene_name", "ensembl_gene_id"),
                     filters = c("biotype"), 
                     values = list("protein_coding"),
                     mart = mart)

    #remove transcripts placed on patches
    genes_info <- genes_info[grep('PATCH',genes_info$chromosome_name,invert=T),]
    #remove HSCHR - alleles
    genes_info <- genes_info[grep('HSCHR',genes_info$chromosome_name,invert=T),]
    
    # after that some genes have several records, i.e. TAP2:
    #6	32789610	32806557	TAP2	ENSG00000204267
    #6	32781544	32806599	TAP2	ENSG00000250264
    
    # sorting by ensembl ID, and using the first one
    genes_info <- genes_info[order(genes_info$external_gene_name,genes_info$ensembl_gene_id),]
    
    genes_info <- genes_info[!duplicated(genes_info$external_gene_name),]
    
    
    genes_info <- genes_info[order(genes_info$chromosome_name,genes_info$start_position),]
    
    write.table(genes_info,
                "protein_coding_genes.bed",
                sep = "\t", quote = F, row.names = F, col.names = F)
}

# input: genes.csv - one column csv file, ensembl_gene_id column
# output: genes.bed, has to be sorted with bedtools
gene_vector2bed <- function(genes.csv){
    mart <- init_mart_human()
    genes.bed <- str_replace(genes.csv, "csv", "bed")
    genes <- read_csv(genes.csv)
    gene_coordinates <- get_gene_coordinates(genes$ensembl_gene_id, mart)
    write_tsv(gene_coordinates, genes.bed, col_names = F)
}

# start and end of the gene, all exons
# input = vector of genes, either ENSEMBL_IDS or external names = disease_panel.list.txt, no header
# output = bed file with coordinates = disease_panel.list.bed
# output is not sorted please sort with bedtools or bash sort
get_gene_coordinates <- function(v_ensembl_gene_ids, mart){
    # test:  
    # https://raw.githubusercontent.com/naumenko-sa/cre/master/data/lupus.csv
    # gene_list_csv <- "lupus.csv"
    # gene_list_csv <- "immunodeficiency.csv"
    
    #guess gene id type
    #gene_ids <- read.csv(gene_list_csv, stringsAsFactors = F)
    
    #assumming that the first column is has ENSEMBL_GENE_IDs
    #agene <- gene_ids[1,1]
    
    #if (grepl("ENSG", agene, fixed=T)){
    #    filter = "ensembl_gene_id"
    #}else{
    #    filter = "external_gene_name"
    #}
        
    genes <- getBM(
        attributes = c("ensembl_gene_id", "chromosome_name", 
                       "start_position", "end_position",
                       "external_gene_name"),
        filters = c("ensembl_gene_id"),
        values = v_ensembl_gene_ids,
        mart = mart)
    
    genes_bed <- as_tibble(genes[c(2:4)])
    
    genes_bed <- genes_bed %>% 
                    filter(str_detect(chromosome_name, "PATCH", negate = T)) %>% 
                    filter(str_detect(chromosome_name, "HSCHR", negate = T))
                                                 
    # sort with bedtools
    #genes_bed <- genes_bed[order(as.numeric(as.character(genes_bed$chromosome_name)),
    #                             genes_bed$start_position),]

    # for the assignment!!! = reference hg19 not grch37
    # genes_bed$chromosome_name <- paste0("chr",genes_bed$chromosome_name)
    return(genes_bed)
}

# reads HPO.tsv file pulled from Phenotips and generates a file with gene coordinates
phenotips_hpo2gene_coordinates <- function(phenotips_hpo.tsv){
    #phenotips_hpo.tsv <- "1153_CH0769_HPO.tsv"
    hpo_genes <- read_tsv(phenotips_hpo.tsv)
    hpo_genes.missing_ensg <- filter(hpo_genes, str_detect(`Gene ID`, "ENSG", negate = T))
    
    cat("Genes:\n")
    print(sort(hpo_genes.missing_ensg$`Gene ID`))
    cat("are missing ensembl_gene_id, please add corresponding intervals manually to the bed file\n")
    cat("Some genes are listed in ~/cre/data/missing_genes_grch37.bed\n")
    missing_genes <- read_tsv("~/cre/data/missing_genes_grch37.bed")
    options(tibble.print_max = Inf)
    print(arrange(missing_genes, gene))
    
    hpo_genes <- hpo_genes %>% filter(str_detect(`Gene ID`, "ENSG")) 
    v_ensembl_gene_ids <- hpo_genes$`Gene ID`
    
    mart <- init_mart_human_grch37()
    hpo_genes_bed <- get_gene_coordinates(v_ensembl_gene_ids, mart)
      
    output_file_name <- gsub(".tsv",".unsorted.bed", phenotips_hpo.tsv)
    
    cat(paste0("Writing gene intervals including UTRs to", output_file_name,"\n"))
    cat("Please sort and merge bed file with bedtools before using in HPC pipelines\n")
    #bedtools does not like bed headers
    write_tsv(hpo_genes_bed, output_file_name, col_names = F)
}
    

# sometimes people want ccds genes, then use with_ccds
# https://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi
# some genes don't have CCDS while they are coding, i.e. B4GALT1, ISPD, LARGE
get_ccds_genes_chr <- function(chromosome, mart){
    #test
    chromosome <- "X"
    ccds_genes_chr <- getBM(attributes = c("ensembl_gene_id", "ccds"),
                            filters = c("chromosome_name","with_ens_hs_translation"),                     
                            values = list(chromosome, T),
                     mart = mart)
    
    return(ccds_genes_chr)
    
    # ccds = ccds id
}

#LSP1 corresponds to two genes and has exons on chr11 and chr13 - a bug - fixed in GRCH38
#CKS1B: chr1 and chr5
tutorial_lsp1_gene <- function(){
    #mart <- init_mart_human()
    mart_grch38 <- tutorial_init_mart_human_grch38()
    lsp1_bug <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id",
                                "transcript_count", "ensembl_exon_id",
                                "chromosome_name", "exon_chrom_start",
                                "exon_chrom_end", "genomic_coding_start",
                                "genomic_coding_end", "external_gene_name"),
                  filters = c("external_gene_name"),
                  values = list("LSP1"),
                  mart=mart_grch38)
}

# use chromosomes because of biomart webservice's timeout
# tutorial: 
# - what is exon_chrom_start, or genomic_coding_start?
# - attributes = listAttributes(mart, what = c("name", "description", "fullDescription"))
# - noncoding exons
get_exon_coordinates_chr <- function(chromosome, mart){
    # test:
    # chromosome = 21
    genes_for_chr <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "transcript_count",
                                          "ensembl_exon_id", "chromosome_name", "exon_chrom_start",
                                          "exon_chrom_end", "genomic_coding_start", "genomic_coding_end",
                                          "external_gene_name"),
                        filters = c("chromosome_name"),
                        values = list(chromosome),
                        mart = mart) %>% 
                     as_tibble() %>% 
                     mutate_at(vars(chromosome_name), as.character)
    
    return(genes_for_chr)
}

# print genomic_coding coordinates and exclude UTRs
# 18710 genes
# result it is neither sorted not merged!
# in bash:
# cat coding.exons.bed | sort -k1,1 -k2,2n > coding.exons.sorted.bed
# bedtools merge -i coding.exons.bed.sorted -c 4 -o distinct > coding.exons.merged.bed
get_exon_coordinates <- function(mart){
    
    exon_coordinates <- NULL
    # MT is problematic in ccds
    # ~5 minutes for all chromosomes
    for (chr in c(seq(1, 22), "X", "Y")){
        print(chr)
        buffer <- get_exon_coordinates_chr(chr, mart)
        exon_coordinates <- bind_rows(exon_coordinates, buffer)
    }
    
    # remove noncoding exons
    # select is masked by biomart
    # remove 1bp exons
    exon_coordinates <- exon_coordinates %>%
        drop_na(genomic_coding_start, genomic_coding_end) %>% 
        dplyr::select(chromosome_name, genomic_coding_start, genomic_coding_end, external_gene_name) %>% 
        filter(genomic_coding_start != genomic_coding_end) %>% 
        write_delim("coding_exons.unsorted.bed", col_names = F)
    
    #library("bedr")
    #exon_coordinates.bed.merged = bedr(input = list(i=exon_coordinates.bed.sorted),
    #                                   method="merge",engine="bedtools",check.chr = F, check.zero.based = F,
    #                                   params = "-c 4 -o distinct")
    
    #write.table(exon_coordinates.bed.merged,"coding.exons.bed",sep="\t",quote=F,row.names=F,col.names=F)
    
    #IRanges: within a chromosome!
    #exon_coordinates.range = IRanges(start = genomic_coding_start,
    #                                 end = genomic_coding_end,
    #                                 names = external_gene_name)
    
    #exon_coordinates.range = reduce(exon_coordinates.range)
    #exon_coordinates.bed = as.data.frame(exon_coordinates.range)
    
    #for individual genes
    #for (gene in sort(unique(exon_coordinates.bed$external_gene_name)))
    #{
    #    print(gene)
    #    gene.table = exon_coordinates.bed[exon_coordinates.bed$external_gene_name == gene,]
    #    gene.bed = subset(gene.table, 
    #                      select=c('chromosome_name','genomic_coding_start',
    #                               'genomic_coding_end','ensembl_exon_id'))
    #    gene.bed = gene.bed[order(gene.bed$genomic_coding_start),]
        
        #have to sort and merge this with bedtools
    #    write.table(gene.bed,paste0(gene,".unsorted.bed"),sep='\t',quote=F,row.names=F,col.names=F)
    #}
}

tutorial_get_sequence <- function(){
    #seq = getSequence(id="ENST00000357033",
    #                  type="ensembl_transcript_id",
    #                  seqType = "coding", mart=mart)
    seq <- getSequence(id = "ENSG00000172062",
                       type = "ensembl_gene_id",
                       seqType = "gene_exon_intron", 
                       mart = mart)
    
    
    write(">SMN1", "SMN1.fasta")
    write(seq$gene_exon_intron, "SMN1.fasta", append = T)
    
    peptide <- getSequence(id = "ENSG00000172062",
                          type = c("ensembl_gene_id"),
                          seqType = "peptide", 
                          mart = mart)
}

bedtools_sort_and_merge_example <- function(){
  # installation:
  # install bedtools and bedops and put them to your PATH
  # http://bedtools.readthedocs.io/en/latest/content/installation.html
  # https://github.com/bedops/bedops/releases
  # documentation:
  # https://cran.r-project.org/web/packages/bedr/bedr.pdf
  
  index <- get.example.regions()
  a = index[[1]]
  a.sorted = bedr(engine="bedtools",input = list(i=a), method="sort", params="")
  a.merged = bedr(engine="bedtools",input = list(i=a.sorted), method="merge", params="")
  
}

# input file with gene_name column = panel.genes.csv
# output file = panel.csv
get_ensembl_gene_ids <- function(panel.genes.csv){
    mart <- init_mart_human()
    genes <- read_csv(panel.genes.csv)
    
    ensembl_ids <- get_ensembl_gene_ids2(genes$gene_name, mart)
    
    panel.csv <- str_replace(panel.genes.csv, ".genes", "")
    write_excel_csv(ensembl_ids, panel.csv)
}

get_ensembl_gene_ids2 <- function(v_gene_names, mart, keep_duplicates = F){
    # test:
    # v_gene_names_file = "omim.gene.list"
    ensembl_genes <- as_tibble(getBM(attributes = c("ensembl_gene_id", "external_gene_name",
                                                    "chromosome_name"),
              filters = c("external_gene_name"),
              values = v_gene_names,
              mart = mart))
    
    # sometimes genes on alleles have lesser ENSG than main assembly (CDSN in grhc37)
    ensembl_genes <- ensembl_genes %>% 
        filter(str_detect(chromosome_name, "PATCH", negate = T)) %>% 
        filter(str_detect(chromosome_name, "HSCHR", negate = T))
    
    ensembl_genes$chromosome_name <- NULL
    
    genes_not_found <- v_gene_names[!(v_gene_names %in% ensembl_genes$external_gene_name)]
    if (length(genes_not_found) > 0 ){
        cat("ensembl_gene_id not found for: \n")
        print(genes_not_found)
        cat("trying to fetch them from grch38 \n")
        mart38 <- init_mart_human(host = "www.ensembl.org")
        ensembl_genes38 <- as_tibble(getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                                         filters = c("external_gene_name"),
                                         values = genes_not_found,
                                         mart = mart38))
        still_not_found <- genes_not_found[!(genes_not_found %in% ensembl_genes38$external_gene_name)]
        if (length(still_not_found) > 0){
            cat("Still can't pull genes: \n")
            print(still_not_found)
        }
        ensembl_genes <- bind_rows(ensembl_genes, ensembl_genes38)
    }
    
    ensembl_genes <- ensembl_genes %>% 
        arrange(external_gene_name, ensembl_gene_id)
    
    if (keep_duplicates == T){
        #keeping all ensembl IDs for manual curation
        ensembl_genes <- ensembl_genes %>% 
            group_by(external_gene_name) %>% 
            summarise(ensembl_gene_id = paste(ensembl_gene_id, collapse =","))
    }else{
        # if there are multiple ENS_ID for a gene name, we keep the smallest ENS_ID
        # which is not alsways the case:
        # C4A: ENSG00000244731 not ENSG00000244207
        ensembl_genes <- ensembl_genes %>% distinct(external_gene_name, .keep_all = T)
    }
    
    return(ensembl_genes)
}

get_external_gene_names <- function(gene_list_file){
    gene_list_file <- "omim.gene.list"
    genes <- read.table(gene_list_file,stringsAsFactors=F)
    genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                filters = c("ensembl_gene_id"),
                values = genes, mart = mart)
    write.table(genes,paste0(gene_list_file, ".w_names"), sep = "\t", quote = F, 
                row.names = F, col.names = F)
}

# coordinates of the exon starts and ends 
# some exons of the canonical isoform are non coding - they have NA in genomic_coding_start,
# they are excluded, for example in ACTA1
# -001 is NOT always the canonical isoform
# takes the longest cds from gencode_basic transcripts
# PLEC gene has two ensembl identifiers:
# ENSG00000178209 - we need this one, which is protein_coding
# ENSG00000261109 
# this is a slow function, use it for individual genes/gene panels, when every exon is needed
# for bulk coverage analysis use the function above
get_exon_coordinates_for_canonical_isoform <- function(ensembl_gene_id, mart){
    # gene_name='HNRNPDL'
    # gene_name = 'PLEC'
    # PATCH gene
    # gene_name = 'ABBA01057584.1'
    # gene_name = 'SEPN1'
    
    # a problematic gene, ENSEMBL returns it on HSCHR6_MHC_COX if you set biotype filter = protein_coding, it will return HSCHR6
    # gene_name='VARS2'
    
    # has NA in CDS_length, does not have genomic coding start and end
    # gene_name="RMRP" 
    # gene_name="SDHAF2"
    # gene_name="MAX"
    # gene_name = "SLC7A7"
    
    # test
    # ensembl_gene_id <- "ENSG00000000419"
    cat("Processing gene:", ensembl_gene_id, ",")
    genes_info <- as_tibble(getBM(attributes = c("chromosome_name", "external_gene_name", "ensembl_transcript_id",
                                       "cds_length", "ensembl_gene_id"),
                        filters = c("ensembl_gene_id", "transcript_gencode_basic"), 
                        values = list(ensembl_gene_id = ensembl_gene_id,
                                      transcript_gencode_basic = T),
                        mart = mart))
                     
    #remove transcripts located on patches and alleles (HSCHR)
    genes_info <- genes_info %>% 
                    filter(!grepl("PATCH", chromosome_name)) %>% 
                    filter(!grepl("HSCHR", chromosome_name))
    
    # select canonical transcript print out the single trancript
    genes_info <- arrange(genes_info, desc(cds_length)) 
    canonical_transcript <- genes_info$ensembl_transcript_id[1]
        
    cat(" Canonical transcript:", canonical_transcript, "\n")
    
    if (is.na(canonical_transcript)){ cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n
No transcript for ", ensembl_gene_id, "\n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n")}
    
    genes_info <- as_tibble(getBM(attributes = c("chromosome_name", 
                                       "genomic_coding_start", "genomic_coding_end",
                                       "ensembl_exon_id", "external_gene_name", "ensembl_gene_id",
                                       "start_position", "end_position",
                                       "exon_chrom_start", "exon_chrom_end", 
                                       "ensembl_transcript_id"),
                filters = c("ensembl_transcript_id"),
                values = c(canonical_transcript), mart = mart)) %>% 
        mutate(chromosome_name = as.character(chromosome_name)) %>%
        filter(!is.na(genomic_coding_start) & !is.na(genomic_coding_end)) %>%
        select(chromosome_name, genomic_coding_start, genomic_coding_end, external_gene_name)
    return(genes_info)
}

# get exon coordinates for canonical isoform for genes in a list
# input: genes.csv - 1 column csv with ensembl_gene_id
# output: genes.bed
get_exon_coordinates2 <- function(genes_csv){
    mart <- init_mart_human()
    genes <- read_csv(genes_csv)
    
    exons <- get_exon_coordinates_for_canonical_isoform(genes$ensembl_gene_id[1], mart)
    
    for(gene in tail(genes$ensembl_gene_id, -1)){
        exons_buf <- get_exon_coordinates_for_canonical_isoform(gene, mart)
        exons <- bind_rows(exons, exons_buf)
    }
    
    genes_bed <- str_replace(genes_csv, "csv", "bed")
    write_tsv(exons, genes_bed, col_names = F)
}

#for (gene in c("SETX","PNKP","AP3B2","GUF1"))
#{
#    get_exon_coordinates_for_canonical_isoform(gene,mart)
#}

#exon coordinates given ENS ids
get_omim_orphanet_exon_coordinates <- function(){ 
    omim_orphanet_ens_ids <- read.table("omim.orphanet.v2.ENS")
  
    omim_exons <- getBM(attributes = c('ensembl_gene_id','ensembl_transcript_id','transcript_count',
                                       'ensembl_exon_id','chromosome_name','exon_chrom_start','exon_chrom_end',
                                       'genomic_coding_start', 'genomic_coding_end'),
                    filters = c('ensembl_gene_id'),
                    values = omim_orphanet_ens_ids, mart = grch37)
  
    omim_exons <- na.omit(omim_exons)
  
    omim_exons.short <- as.data.frame(unique(omim_exons[c("ensembl_gene_id","chromosome_name")]))
    omim_exons.short.table <- as.data.frame.table(table(omim_exons.short$chromosome_name))
    sum(omim_exons.short.table[grep("HG|HS",omim_exons.short.table$Var1),]$Freq)
  
    omim_exons_chr <- omim_exons[grep("HG|HS",omim_exons$chromosome_name,invert=T),]
    
    omim_exons.grch38 <- getBM(
    attributes <- c('ensembl_gene_id','ensembl_transcript_id','transcript_count','ensembl_exon_id',
                 'chromosome_name','exon_chrom_start','exon_chrom_end'),
    filters <- c('ensembl_gene_id'),
    values <- omim_orphanet_ens_ids,mart=grch38)
    omim_exons.grch38.short <- as.data.frame(unique(omim_exons.grch38[c("ensembl_gene_id","chromosome_name")]))
    omim_exons.grch38.short.table <- as.data.frame.table(table(omim_exons.grch38.short$chromosome_name))
    sum(omim_exons.grch38.short.table[grep("HG|HS", omim_exons.short.table$Var1),]$Freq)

    write.table(omim_exons_chr, "omim.exons", quote=F, row.names=F, col.names=F)
    omim_exons.bed=subset(omim_exons_chr, select = c("chromosome_name","genomic_coding_start","genomic_coding_end"))
    write.table(omim_exons.bed,"omim.exons.notsorted.bed",sep="\t",quote=F,row.names=F,col.names=F)
}

# read counts for all genes, filter protein coding genes only
filter_protein_coding_genes <- function(counts.csv){
    counts <- read_csv(counts.csv)
    protein_coding_genes <- read_csv("~/cre/data/protein_coding_genes.csv")
    counts <- inner_join(counts, protein_coding_genes, by = c("Ensembl_gene_id" = "ENS_GENE_ID"))
    output.csv <- str_replace(counts.csv, "csv", "coding.csv")
    write_excel_csv(counts, output.csv)
}

###############################################################################
args <- commandArgs(trailingOnly = T)
if (length(args) == 0 || args[1] == "--help"){
    cat("Usage: Rscript function_name function_args\n")
    cat("Available functions:\n")
    cat("phenotips_hpo2gene_coordinates phenotips_hpo.tsv\n")
    cat("gene_vector2bed genes.csv\n")
    cat("get_exon_coordinates2 genes.csv\n")
    cat("get_ensembl_gene_ids panel.genes.csv\n")
    cat("filter_protein_coding_genes counts.csv\n")
}else{
    cat(paste0("Running function: ", args[1],"\n"))
    init()
    fcn <- get(args[1])
    fcn(tail(args,-1))
}
###############################################################################