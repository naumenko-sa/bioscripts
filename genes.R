# biomart wrappers to get gene,transcript,exons annotations from ENSEMBL

init = function()
{
    library("biomaRt")  
    mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", 
                   path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
    datasets=listDatasets(mart)
  
    mart = useDataset(mart,dataset="hsapiens_gene_ensembl")
  
    attributes=listAttributes(mart)
    filters=listFilters(mart)
  
    chromosomes = getBM(attributes=c('chromosome_name'),mart=mart)
    
    return(mart)
    
    #grch38 = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
    #datasets = listDatasets(grch38)
    #grch38 = useDataset(grch38,dataset="hsapiens_gene_ensembl")
}

get_protein_coding_genes = function(mart)
{
    setwd("~/Desktop/reference_tables/")
    protein_coding_genes = getBM(attributes=c('ensembl_gene_id','external_gene_name'),
                                 filters = c('biotype'),
                                 values='protein_coding',
                                 mart=mart)
    colnames(protein_coding_genes)[2] = 'gene_name'
    #genes might be not unique - polymorphic regions like NCR3 gene
    write.table(unique(sort(protein_coding_genes[,2])),file="protein_coding_genes.list",quote=F,row.names=F, col.names=F,sep="\t")
}

get_gene_descriptions = function()
{
    #name_1006 is GO_term
    ensembl_w_description = getBM(attributes=c('ensembl_gene_id','external_gene_name','description','name_1006'),mart=grch37)
    write.table(ensembl_w_description,file="ensembl_w_description.txt1",quote=F,row.names=F,sep="\t")
}

get_refseq_transcript_ids = function()
{
    refseq_transcripts = getBM(attributes=c('ensembl_transcript_id','refseq_mrna'),mart=grch37)
    refseq_transcripts = getBM(attributes=c('ensembl_transcript_id','external_transcript_name'),mart=grch37)
    
    write.table(refseq_transcripts[refseq_transcripts$refseq_mrna!='',],file="ensembl_refseq.txt",quote=F,row.names=F,sep="\t")
}

#use chromosomes because of biomart webservice timeout
#sometimes people want ccds genes, then use with_ccds
#some genes don't have CCDS while they are coding, i.e. B4GALT1, ISPD, LARGE
get_exon_coordinates_chr = function(chromosome)
{
    #ccds_genes = getBM(attributes=c('ensembl_gene_id'),
    #                 filters=c('with_ens_hs_translation','chromosome_name'),
    #                 values=list(T,chromosome),
    #                 mart=grch37)
    
    #getBM(
    #  attributes=c('ensembl_gene_id','ensembl_transcript_id','transcript_count','ensembl_exon_id',
    #    'chromosome_name','exon_chrom_start','exon_chrom_end','genomic_coding_start','genomic_coding_end',
    #    'external_gene_name'),
    #  filters=c('ensembl_gene_id'),
    #  values=list(ccds_genes),
    #  mart=grch37)
  
    getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','transcript_count','ensembl_exon_id',
          'chromosome_name','exon_chrom_start','exon_chrom_end','genomic_coding_start','genomic_coding_end',
          'external_gene_name'),
          filters = c('chromosome_name'),
          values = list(chromosome),
          mart=grch37)
}

#LSP1 has exons on chr11 and chr13 - a bug to report
#CKS1B: chr1 and chr5
test_lsp1_gene = function()
{
  library("biomaRt")  
  grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", 
                   path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
  datasets=listDatasets(grch37)
  
  grch37 = useDataset(grch37,dataset="hsapiens_gene_ensembl")
  
  attributes=listAttributes(grch37)
  filters=listFilters(grch37)
  
  lsp1_bug=getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','transcript_count','ensembl_exon_id',
                     'chromosome_name','exon_chrom_start','exon_chrom_end','genomic_coding_start','genomic_coding_end',
                     'external_gene_name'),
        filters = c('external_gene_name'),
        values = list("LSP1"),
        mart=grch37)
}

#print genomic_coding to exclude UTRs
get_exon_coordinates = function()
{
    #18710 genes
    exon_coordinates=get_exon_coordinates_chr(1)
    #MT is problematic in ccds
    for (chr in c(seq(2,22),"X","Y"))
    {
        buffer = get_exon_coordinates_chr(chr)
        exon_coordinates=rbind(buffer,exon_coordinates)
    }
    exon_coordinates=na.omit(exon_coordinates)
    write.table(exon_coordinates,"ccds.coding.exons",quote=F,row.names=F,col.names=F)
    write.table(unique(exon_coordinates$ensembl_gene_id),"ccds.coding.genes.ENS",quote=F,row.names=F,col.names=F)
    exon_coordinates.bed=subset(exon_coordinates,select=c("chromosome_name","genomic_coding_start","genomic_coding_end","external_gene_name"))
    write.table(exon_coordinates.bed,"ccds.coding.exons.notsorted.bed",sep="\t",quote=F,row.names=F,col.names=F)
    
    ccds_genes = getBM(attributes=c('ensembl_gene_id'),mart=grch37)
}

# coordinates of the gene start and gene end (all exons)
get_gene_coordinate = function(gene_list_file)
{
    genes = read.table(gene_list_file,stringsAsFactors=F)
    genes=getBM(
      attributes=c('ensembl_gene_id','chromosome_name','start_position','end_position','external_gene_name'),
      filters=c('external_gene_name'),
      values=genes,mart=grch37)
    write.table(genes[c(2:5)],paste0(gene_list_file,".bed"),sep="\t",quote=F,row.names=F,col.names=F)
}

# coordinates of the exon starts and ends 
# some exons of the canonical isoform are non coding - they have NA in genomic_coding_start,
# they are excluded, for example in ACTA1
# -001 is NOT always the canonical isoform
# takes the longest cds from gencode_basic transcripts
# PLEC gene has two ensembl identifiers:
# ENSG00000178209 - we need this one, which is protein_coding
# ENSG00000261109 
get_exon_coordinates_for_canonical_isoform = function(gene_name,mart)
{
    #gene_name='HNRNPDL'
    #gene_name = 'PLEC'
    print(gene_name)
    genes_info=getBM(attributes=c('chromosome_name','external_gene_name','ensembl_transcript_id','cds_length','ensembl_gene_id'),
                     filters=c('external_gene_name','transcript_gencode_basic','biotype'), 
                     values=list(external_gene_name=gene_name,transcript_gencode_basic=T,'protein_coding'),mart=mart)
    
    #remove transcripts placed on patches
    genes_info = genes_info[grep('PATCH',genes_info$chromosome_name,invert=T),]
    
    if (nrow(genes_info>1)){
        genes_info = genes_info[order(-genes_info$cds_length),]  
    }
    canonical_transcript = genes_info$ensembl_transcript_id[1]
    
    genes_info=getBM(attributes=c('chromosome_name','genomic_coding_start','genomic_coding_end','ensembl_exon_id',
                                  'external_gene_name','ensembl_gene_id',
                                  'start_position','end_position',
                                  'exon_chrom_start','exon_chrom_end','ensembl_transcript_id'),
                    filters=c('ensembl_transcript_id'), values=c(canonical_transcript),mart=mart)
    
    genes_info = na.omit(genes_info)
    
    #bedtools does not like colnames
    write.table(genes_info[c(1:4)],paste0(gene_name,".bed"),sep="\t",quote=F,row.names=F,col.names=F)
    write.table(genes_info,paste0(gene_name,".extended.bed"),sep="\t",quote=F,row.names=F,col.names=T)
    print(paste(gene_name,genes_info[1,4]))
}

get_exon_coordinates_for_muscular_genes = function()
{
    mart=init()
    setwd("~/Desktop/project_RNAseq_diagnostics/gene_panels/")
    muscular_gene_panels = read.csv("/home/sergey/Desktop/project_RNAseq_diagnostics/gene_panels/muscular_gene_panels.csv", stringsAsFactors=F)
    
    for(gene in unique(sort(muscular_gene_panels$Gene)))
    {
        get_exon_coordinates_for_canonical_isoform(gene,mart)
    }
}

#exon coordinates given ENS ids
get_omim_orphanet_exon_coordinates = function()
{ 
  omim_orphanet_ens_ids = read.table("omim.orphanet.v2.ENS")
  
  omim_exons=getBM(
      attributes=c('ensembl_gene_id','ensembl_transcript_id','transcript_count','ensembl_exon_id',
      'chromosome_name','exon_chrom_start','exon_chrom_end','genomic_coding_start','genomic_coding_end'),
      filters=c('ensembl_gene_id'),
      values=omim_orphanet_ens_ids,mart=grch37)
  
  omim_exons=na.omit(omim_exons)
  
  omim_exons.short = as.data.frame(unique(omim_exons[c("ensembl_gene_id","chromosome_name")]))
  omim_exons.short.table = as.data.frame.table(table(omim_exons.short$chromosome_name))
  sum(omim_exons.short.table[grep("HG|HS",omim_exons.short.table$Var1),]$Freq)
  
  omim_exons_chr=omim_exons[grep("HG|HS",omim_exons$chromosome_name,invert=T),]
    
  omim_exons.grch38=getBM(
    attributes=c('ensembl_gene_id','ensembl_transcript_id','transcript_count','ensembl_exon_id',
                 'chromosome_name','exon_chrom_start','exon_chrom_end'),
    filters=c('ensembl_gene_id'),
    values=omim_orphanet_ens_ids,mart=grch38)
  omim_exons.grch38.short = as.data.frame(unique(omim_exons.grch38[c("ensembl_gene_id","chromosome_name")]))
  omim_exons.grch38.short.table = as.data.frame.table(table(omim_exons.grch38.short$chromosome_name))
  sum(omim_exons.grch38.short.table[grep("HG|HS",omim_exons.short.table$Var1),]$Freq)

  write.table(omim_exons_chr,"omim.exons",quote=F,row.names=F,col.names=F)
  omim_exons.bed=subset(omim_exons_chr,select=c("chromosome_name","genomic_coding_start","genomic_coding_end"))
  write.table(omim_exons.bed,"omim.exons.notsorted.bed",sep="\t",quote=F,row.names=F,col.names=F)
  
}

setwd("~/Desktop/tools/MendelianRNA-seq/data/")
mart=init()
get_exon_coordinates_for_canonical_isoform("DMD",mart)

get_gene_coordinate("kidney.glomerular.genes")

get_exon_coordinates()
get_omim_orphanet_exon_coordinates()

setwd("~/Desktop/reference_tables/")
get_gene_coordinate("protein_coding_genes")

#better to use ENS ids from OMIM/Orphanet text files
#ccds_omim_genes = getBM(attributes=c('ensembl_gene_id','mim_gene_accession','mim_morbid_accession'),
#        filters=c('with_ccds','with_mim_gene','with_mim_morbid'),
#        values=list(T,T,T),
#        mart=grch37)