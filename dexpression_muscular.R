library(edgeR)
setwd("~/Desktop/project_muscular/")

all_counts=read.delim("annotated_combined.counts")
#x=all_counts[c("blood1","blood2","fibroblast5","muscle1","muscle2","myotubes5")]
x=all_counts[c("muscle1","muscle2")]
group=factor(c(1,2))
  
output="muscle1.muscle2"
  
y=DGEList(counts=x,group=group)
  
  #filtration
  #keep=rowSums(cpm(y)>2)==35
keep=rowSums(cpm(y)>1) ==2
y=y[keep,,keep.lib.sizes=FALSE]
  
  #normalized counts
nc=cpm(y,normalized.lib.sizes=F)
write.table(nc,paste0(output,".normalized_counts.txt"),col.names=NA)
  
  png(paste0(output,".mds.raw.png"))
  plotMDS(y,las=1,main=paste0(output," - MDS plot for raw data"))
  dev.off()
  
  y=calcNormFactors(y)
  
  png(paste0(output,".mds.norm.png"))
  plotMDS(y,las=1,main=paste0(output," - MDS plot for normalized data"))
  dev.off()
  
  write.table(y$samples,paste0(output,".norm.factors.txt"),col.names=NA)
  
  design=model.matrix(~group)
  y=estimateDisp(y,design)
  
  fit=glmFit(y,design)
  lrt=glmLRT(fit)
  #prints top 50 genes with p.value<0.05, check if there are more
  write.table(topTags(lrt,p.value=0.05,n=50),paste0(output,".significant_genes.txt"),col.names=NA)
  write.table(lrt$table,paste0(output,".all_genes.txt"),col.names=NA)
}




group=factor(c(1,1,1,2,2,2))
output="test1"
y=DGEList(counts=x,group=group)
y.cpm=cpm(y)

#filtration
#keep=rowSums(cpm(y)>2)==35
keep=rowSums(cpm(y)>1) >=4
y=y[keep,,keep.lib.sizes=FALSE]

#mdsPlot
colors=substr(colnames(y$counts),1,3)
colors=gsub("X12","red",colors)
colors=gsub("X1T","green",colors)
colors=gsub("X6T","blue",colors)

plotMDS(y,las=1,col=colors,main="MDS plot muscular dataset",xlim=c(-0.4,0.4),ylim=c(-0.4,0.4))
plotMDS(y,las=1,main="MDS plot for muscular dataset")


nc=cpm(y,normalized.lib.sizes=F)
write.table(nc,"filtered.normalized_counts.txt",col.names=NA)

plotMDS(nc,las=1,main="MDS plot for CPM") 

y=calcNormFactors(y)
write.table(y$samples,"filtered.norm.factors.txt",col.names=NA)    
plotMDS(y,las=1,col=colors,main="MDS plot for normalized data",ylim=c(0.4,-0.4))
plotMDS(y,las=1,main="MDS plot for normalized data")


#FXN
fxn=
