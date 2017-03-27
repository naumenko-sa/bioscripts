#basic function for dexpression

#merges individual count files into 1
merge_counts = function()
{
  
  samples = read.table("samples.txt", quote="\"", stringsAsFactors=F)
  samples = samples[,1]
  
  for (sample in samples)
  {
    table_name = sample
    file_name = paste0(sample,".counts")
    assign(table_name,read.delim(file_name, header=F,row.names=1,stringsAsFactors=F))
    table_itself = get(table_name)
    colnames(table_itself) = c(sample)
    assign(table_name,table_itself)
  }
  
  samples.data = get(samples[1])
  
  for (sample in tail(samples,-1))
  {
    #test
    #table_name = samples[2]
    table_name = sample
    table_itself = get(table_name)
    samples.data = merge(samples.data,get(table_name),by.x="row.names",by.y="row.names")
    row.names(samples.data)=samples.data$Row.names
    samples.data$Row.names = NULL
  }
  
  write.table(samples.data,"all_counts.txt",quote=F)
}