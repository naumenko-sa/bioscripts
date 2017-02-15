# Loads a file with counts from feature_counts
# loads gene lengths
# loads gene names
# calculates RPKMs
# returns ENS_ID, rpkm, Gene_name
load_rpkm_counts = function(filename)
{
    #test:
    #filename="/home/sergey/Desktop/project_muscular/Fibroblast8/fibroblast8.rpkm"
    library(edgeR)   
    ensembl_w_description = read.delim2("~/Desktop/reference_tables/ensembl_w_description.txt", row.names=1, stringsAsFactors=F)
    #first line in the file is a comment
    counts = read.delim(filename, stringsAsFactors=F, row.names=1,skip=1)
    counts$Chr=NULL
    counts$Start=NULL
    counts$End=NULL
    counts$Strand=NULL
    
    
    Gene_lengths = counts$Length
    
    counts$Length=NULL
    
    counts = rpkm(counts,Gene_lengths)
    
    counts = merge(counts,ensembl_w_description,by.x="row.names",by.y="row.names")
    row.names(counts)=counts$Row.names
    counts$Row.names=NULL
    counts$Gene_description=NULL
    
    return(counts)
}

# merge two dataframes by row.names and fix the row.names of the resulting df
merge_row_names = function(df1,df2)
{
  merged = merge(df1,df2,by.x='row.names',by.y='row.names')
  row.names(merged) = merged$Row.names
  merged$Row.names = NULL
  return(merged)
}

