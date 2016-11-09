library(edgeR)
setwd("~/Desktop/project_muscular/")

all_counts = read.delim("annotated_combined.counts", row.names=1, stringsAsFactors=FALSE)
norm_counts = cpm(all_counts[,c(1,2,3,4,5,6)])
gene_panel = read.delim("muscle2_genes.txt",stringsAsFactors = F,header=T)

muscular_genes = all_counts[all_counts$symbol %in% unlist(gene_panel),]

work_counts=muscular_genes
attach(work_counts)

work_counts$muscle1_norm = with(work_counts,1000000*muscle1/28641998);
work_counts$muscle2_norm = with(work_counts,1000000*muscle2/60574350);

x=muscular_genes[c("muscle1","muscle2")]
group=factor(c(1,2))
y=DGEList(counts=x,group=group)
y=calcNormFactors(y)

bcv=0.2
#from A.thaliana experiment
dispersion=0.04
  
  design=model.matrix(~group)
  fit=glmFit(y,design,dispersion)
  lrt=glmLRT(fit,coef=2)

write.table(merge(work_counts,lrt$table,by="row.names"),"muscular_panel_counts.txt",col.names=NA,quote=F)  

genes = read.delim("test.txt",stringsAsFactors = F,header=F,sep=",")

