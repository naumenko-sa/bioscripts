# biomart wrappers to get gene information, exon coordinates etc

init = function()
{
    library("biomaRt")  
    grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", 
                   path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
    datasets=listDatasets(grch37)
  
    grch37 = useDataset(grch37,dataset="hsapiens_gene_ensembl")
  
    attributes=listAttributes(grch37)
    filters=listFilters(grch37)
  
    chromosomes = getBM(attributes=c('chromosome_name'),mart=grch37)
    
    grch38 = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
    datasets = listDatasets(grch38)
    grch38 = useDataset(grch38,dataset="hsapiens_gene_ensembl")
}

get_protein_coding_genes = function()
{
    setwd("~/Desktop/reference_tables/")
    protein_coding_genes = getBM(attributes=c('ensembl_gene_id','external_gene_name'),
                                 filters = c('biotype'),
                                 values='protein_coding',
                                 mart=grch37)
    colnames(protein_coding_genes)[2] = 'gene_name'
    write.table(protein_coding_genes,file="protein_coding_genes.txt",quote=F,row.names=F,sep="\t")
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
#some genes don't have CCDS while they are coding, i.e. B4GAT1, ISPD, LARGE1
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
    write.table(unique(exon_coordinates$ensembl_gene_id),"ccds.codins.genes.ENS",quote=F,row.names=F,col.names=F)
    exon_coordinates.bed=subset(exon_coordinates,select=c("chromosome_name","genomic_coding_start","genomic_coding_end","external_gene_name"))
    write.table(exon_coordinates.bed,"ccds.coding.exons.notsorted.bed",sep="\t",quote=F,row.names=F,col.names=F)
    
    ccds_genes = getBM(attributes=c('ensembl_gene_id'),mart=grch37)
}

get_gene_coordinate = function()
{
    #gene_name = "DMD"
    setwd("~/Desktop/project_muscular/reference/")
    muscular_genes = read.table("~/Desktop/project_muscular/reference/muscular_genes.txt", 
                                quote="\"", comment.char="", stringsAsFactors=F)
    genes=getBM(
      attributes=c('ensembl_gene_id','chromosome_name','start_position','end_position','external_gene_name'),
      filters=c('external_gene_name'),
      values=muscular_genes,mart=grch37)
    write.table(genes[c(2:5)],"muscular_genes_coord.bed",sep="\t",quote=F,row.names=F,col.names=F)    
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

setwd("~/Desktop/project_cheo")
init()
get_exon_coordinates()
get_omim_orphanet_exon_coordinates()

#better to use ENS ids from OMIM/Orphanet text files
#ccds_omim_genes = getBM(attributes=c('ensembl_gene_id','mim_gene_accession','mim_morbid_accession'),
#        filters=c('with_ccds','with_mim_gene','with_mim_morbid'),
#        values=list(T,T,T),
#        mart=grch37)