source("~/Desktop/bioscripts/rnaseq.load_rpkm_counts.R")

# calculate RPKM expression of RYR1
setwd("~/Desktop/project_mh/2017-02-01_RYR1_expression/")

samples = unlist(read.table("samples.txt", quote="\"", stringsAsFactors=F))

v=array()
i=1
for(sample in samples)
{
    v[i] = subset(load_rpkm_counts(paste0(sample,".counts.rpkm")),external_gene_name=="RYR1")[[1]]
    i=i+1
}

#GTEX
v[i] = 188.8
samples[i] = "GTEX"

png("naumenko.project_mh.dataset3_rna.RYR1_expression.png",width=1000)
barplot(v,names.arg=samples,main = "Expression of RYR1 gene in MH samples and GTEX, rpkm")
dev.off()