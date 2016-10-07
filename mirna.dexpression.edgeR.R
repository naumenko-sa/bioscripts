#micro RNA dif expression analysis for potato/soybean project

setwd("~/Desktop/project_microRNA_potato/2016-10-07_analysis/")
library(edgeR)

#filter conservative micro RNAs
glymax.conservative = read.table("glymax.conservative", row.names=NULL, stringsAsFactors=FALSE)
glymax.all = read.delim("glymax.tsv", row.names=1, stringsAsFactors=FALSE)
glymax.noncons = setdiff(glymax.all,glymax.conservative)
write.table(glymax.noncons,"glymax.noncons.txt",quote=F,col.names=NA)

soltub.conservative = read.table("soltub.conservative", row.names=NULL, stringsAsFactors=FALSE)
soltub.all = read.delim("soltub.tsv", row.names=1, stringsAsFactors=FALSE)
soltub.noncons = setdiff(soltub.all,soltub.conservative)
write.table(soltub.noncons,"soltub.noncons.txt",quote=F,col.names=NA)

gg=glymax.noncons[,"R4226"]
gs=glymax.noncons[,"R4227"]

ss=soltub.noncons[,"R4227"]
sg=soltub.noncons[,"R4226"]

tt=ss
j=0
for (i in 1:length(tt)) 
  if (tt[i]>0)
    j=j+1
print(j)

plot(soltub.noncons[,"R4226"])

plot(glymax.noncons[,"R4227"])
plot(soltub.noncons[,"R4227"])

g.control=glymax.noncons[c("R4226","R4227")]
s.control=soltub.noncons[c("R4226","R4227")]
x=s.control

#x=all_counts[c("X1T1R1","X1T1R2","X1T1R3","X1T2R1","X1T2R2","X1T2R3")]
#x=all_counts[c("X1T1R1","X1T1R2","X1T1R3","X1T3R1","X1T3R2","X1T3R3")]
x=all_counts[c("R4215", "R4216", "R4217", "R4218")]
x=all_counts[c("R4215", "R4216", "R4217", "R4218","R4219","R4220")]
x=all_counts[c("R4226", "R4227", "R4228", "R4229","R4230","R4231")]


group=factor(c(1,2))
group=factor(c(1,1,2,2))
group=factor(c(1,1,1,1))
group=factor(c(1,1,1,1,1,1))
#group=factor(c(1,1,1,2,2,2))

y=DGEList(counts=x,group=group)

plotMDS(y)
plotMDS(y,cex=1.5,cex.lab=1.5,cex.axis=1.5)

y=calcNormFactors(y)

design=model.matrix(~group)

bcv=0.2
#from A.thaliana experiment
dispersion=0.04

#y=estimateDisp(y,design)

fit=glmFit(y,design,dispersion)
lrt=glmLRT(fit,coef=2)
topTags(lrt,n=262)

top_tags=topTags(lrt,n=262)$table

write.table(top_tags,"glymax.control.txt",quote=F,col.names=NA)
write.table(top_tags,"soltub.control.txt",quote=F,col.names=NA)




o=order(lrt$table$PValue)
normalized_counts=cpm(y)[o[1:262],]
raw_counts=x[o[1:262],]

result=cbind(top_tags,normalized_counts,raw_counts)
filename=paste0(paste(colnames(x),collapse='_'),".txt")
write.table(result,filename)
