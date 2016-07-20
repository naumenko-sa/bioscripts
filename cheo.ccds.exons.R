init = function()
{
  library("biomaRt")  
  grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", 
                   path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
  datasets=listDatasets(grch37)
  
  grch37 = useDataset(grch37,dataset="hsapiens_gene_ensembl")
  
  attributes=listAttributes(grch37)
  filters=listFilters(grch37)
}

get_exon_coordinates_chr = function(chromosome)
{
    ccds_genes = getBM(attributes=c('ensembl_gene_id'),
                     filters=c('with_ccds','chromosome_name'),
                     values=list(T,chromosome),
                     mart=grch37)
    getBM(
      attributes=c('ensembl_gene_id','ensembl_transcript_id','transcript_count','ensembl_exon_id',
        'chromosome_name','exon_chrom_start','exon_chrom_end'),
      filters=c('ensembl_gene_id'),
      values=list(ccds_genes),
      mart=grch37)
}

get_exon_coordinates = function()
{
  #18710 genes
  exon_coordinates=get_exon_coordinates_chr(1)
  for (chr in c(seq(2,22),"X","Y"))
  {
    buffer = get_exon_coordinates_chr(chr)
    exon_coordinates=rbind(buffer,exon_coordinates)
  }
  write.table(exon_coordinates,"ccds.exons",quote=F,row.names=F,col.names=F)
  write.table(unique(exon_coordinates$ensembl_gene_id),"ccds.genes.ENS",quote=F,row.names=F,col.names=F)
  exon_coordinates.bed=subset(exon_coordinates,select=c("chromosome_name","exon_chrom_start","exon_chrom_end"))
  write.table(exon_coordinates.bed,"ccds.exons.notsorted.bed",sep="\t",quote=F,row.names=F,col.names=F)
}

setwd("~/Desktop/project_cheo")
init()
get_exon_coordinates()

#better to use ENS ids from OMIM/Orphanet text files
#ccds_omim_genes = getBM(attributes=c('ensembl_gene_id','mim_gene_accession','mim_morbid_accession'),
#        filters=c('with_ccds','with_mim_gene','with_mim_morbid'),
#        values=list(T,T,T),
#        mart=grch37)