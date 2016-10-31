library(edgeR)
setwd("~/Desktop/project_muscular/")

all_counts=read.delim("annotated_combined.counts")

x=all_counts[c("blood1","blood2","fibroblast5","muscle1","muscle2","myotubes5")]
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
