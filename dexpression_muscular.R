library(edgeR)
setwd("~/Desktop/project_muscular/")

#all_counts = read.delim("annotated_combined.counts", row.names=1, stringsAsFactors=FALSE)
all_counts = read.csv("~/Desktop/project_muscular/all_counts.txt", sep="", stringsAsFactors=FALSE)

symbol = all_counts[,c(1,8)]
#all_counts = all_counts[,c(1,2,3,4,5,6,7,9,10,11,12,13,14,15)]

#c1130_BD_B175 = read.delim("1130-BD-B175.counts", row.names=1,stringsAsFactors=F)

#1258-AC-A79.counts  1275-BK-B225.counts  1388-MJ-M219.counts  2180-CB-C365.counts  6087-PD-B317G.counts  6100-DD-D360.counts


norm_counts = cbind(symbol,cpm(all_counts[,c(2,3,4,5,6,7,9,10,11,12,13,14,15)]))
write.table(norm_counts,"norm_counts_all.txt",col.names=NA,quote=F)  


gene_panel = read.delim("genes_locus.txt",stringsAsFactors = F,header=T)
gene_panel = read.delim("genes_muscular.txt",stringsAsFactors = F,header=T)

muscular_genes = all_counts[all_counts$symbol %in% unlist(gene_panel),]

work_counts=muscular_genes
attach(work_counts)

x=muscular_genes[c("muscle1","X1130.BD.B175")]

x=muscular_genes[c(2,3,4,5,6,7,9,10,11,12,13,14,15)]
group=factor(c(1,2,3,4,5,6,7,8,9,10,11,12,13))

group=factor(c(1,2))

y=DGEList(counts=x,group=group)


y=calcNormFactors(y)

plotMDS(y)
bcv=0.2
#from A.thaliana experiment
dispersion=0.04
  
  design=model.matrix(~group)
  fit=glmFit(y,design,dispersion)
  lrt=glmLRT(fit,coef=2)

write.table(merge(work_counts,lrt$table,by="row.names"),"muscle1_vs_1130_panel.txt",col.names=NA,quote=F)  


