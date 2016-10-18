#micro RNA dif expression analysis for potato/soybean project

setwd("~/Desktop/project_microRNA_potato/2016-10-18-report/")
library(edgeR)

#filter conservative micro RNAs
create_nonconservative_sets = function ()
{
    glymax.conservative = read.table("glymax.conservative", row.names=NULL, stringsAsFactors=FALSE)
    glymax.all = read.delim("glymax.tsv", row.names=1, stringsAsFactors=FALSE)
    glymax.noncons = setdiff(glymax.all,glymax.conservative)
    write.table(glymax.noncons,"glymax.noncons.txt",quote=F,col.names=NA)

    soltub.conservative = read.table("soltub.conservative", row.names=NULL, stringsAsFactors=FALSE)
    soltub.all = read.delim("soltub.tsv", row.names=1, stringsAsFactors=FALSE)
    soltub.noncons = setdiff(soltub.all,soltub.conservative)
    write.table(soltub.noncons,"soltub.noncons.txt",quote=F,col.names=NA)
}

exploratory_analysis = function()
{
  #test
  all.counts = glymax.noncons
  samples = c("R4226","R4227","R4228","R4230")
  samples=colnames(all.counts)
  raw_counts = all.counts[samples]
  
  group=factor(c(1,1,1,1,1,1,1,1,1,1,1,1))
  
  y=DGEList(counts=raw_counts,group=group)
  
  keep = rowSums(cpm(y)>2) >=4
  table(keep)
  
  y = y[keep,,keep.lib.sizes=F]
  
  plotMDS(y)
}

read_nonconservative_sets = function()
{
  glymax.noncons = read.table("glymax.noncons.txt",stringsAsFactors = F)
  soltub.noncons = read.table("soltub.noncons.txt",stringsAsFactors=F)
}

#all.counts = soltub.noncons | glymax.noncons
#background = soltub | glymax
pairwise_comparison = function(sample1,sample2,all.counts,background)
{
  #test
  #sample1="R4227"
  #sample2="R4228"
  #all.counts=soltub.noncons
  #background="soltub"
  
  raw_counts = all.counts[c(sample1,sample2)]
  row.number = nrow(raw_counts)

  group=factor(c(1,2))
  
  y=DGEList(counts=raw_counts,group=group)

  y=calcNormFactors(y)

  design=model.matrix(~group)

  bcv=0.2
  #from A.thaliana experiment
  dispersion=0.04

  #y=estimateDisp(y,design)

  fit=glmFit(y,design,dispersion)
  lrt=glmLRT(fit,coef=2)
  topTags(lrt,n=row.number)

  top_tags=topTags(lrt,n=row.number)

  o=order(lrt$table$PValue)
  normalized_counts=cpm(y)[o[1:row.number],]

  raw_counts=raw_counts[o[1:row.number],]
  
  result=cbind(top_tags,normalized_counts,raw_counts)
  filename=paste0(sample1,"_",sample2,".",background,".result.txt")
  write.table(result,filename,quote=F,col.names=NA)
}

fourfold_comparison = function(samples,all.counts,background)
{
  #test
  samples = c("R4226","R4228","R4227","R4230")
  #all.counts=soltub.noncons
  #background="soltub"
  
  raw_counts = all.counts[samples]
  row.number = nrow(raw_counts)
  
  group=factor(c(1,1,2,2))
  
  y=DGEList(counts=raw_counts,group=group)
  
  y=calcNormFactors(y)
  
  design=model.matrix(~group)
  
  bcv=0.2
  #from A.thaliana experiment
  dispersion=0.04
  
  #y=estimateDisp(y,design)
  
  fit=glmFit(y,design,dispersion)
  lrt=glmLRT(fit,coef=2)
  topTags(lrt,n=row.number)
  
  top_tags=topTags(lrt,n=row.number)
  
  o=order(lrt$table$PValue)
  normalized_counts=cpm(y)[o[1:row.number],]
  
  raw_counts=raw_counts[o[1:row.number],]
  
  result=cbind(top_tags,normalized_counts,raw_counts)
  filename=paste0(paste(samples,collapse="_"),".",background,".result.txt")
  write.table(result,filename,quote=F,col.names=NA)
}

read_nonconservative_sets()
pairwise_comparison("R4226","R4227",soltub.noncons,"soltub")
pairwise_comparison("R4226","R4227",glymax.noncons,"glymax")


pairwise_comparison("R4228","R4227",soltub.noncons,"soltub")
pairwise_comparison("R4226","R4228",soltub.noncons,"soltub")

pairwise_comparison("R4230","R4227",soltub.noncons,"soltub")
pairwise_comparison("R4226","R4230",soltub.noncons,"soltub")

pairwise_comparison("R4227","R4229",glymax.noncons,"glymax")
pairwise_comparison("R4229","R4226",glymax.noncons,"glymax")

pairwise_comparison("R4227","R4231",glymax.noncons,"glymax")
pairwise_comparison("R4231","R4226",glymax.noncons,"glymax")
