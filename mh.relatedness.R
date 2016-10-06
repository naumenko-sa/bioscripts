#calculate relatedness of the families with SNPRelate

#example: https://www.biostars.org/p/83232/
#tutorial: 
#http://www.bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html#format-conversion-from-vcf-files

library(gdsfmt)
library(SNPRelate)
setwd("~/Desktop/project_mh/2016-10-06_families_check/")
#biallelic
snpgdsVCF2GDS("dataset1.vcf", "dataset1.gds")
snpgdsSummary("dataset1.gds")
genofile = snpgdsOpen("dataset1.gds")

#dendogram

dissMatrix  =  snpgdsDiss(genofile , sample.id=NULL, snp.id=NULL, autosome.only=TRUE,remove.monosnp=TRUE, 
                          maf=NaN, missing.rate=NaN, num.thread=2, verbose=TRUE)

snpHCluster =  snpgdsHCluster(dissMatrix, sample.id=NULL, need.mat=TRUE, hang=0.25)

snpgdsClose(genofile)

set.seed(100)
cutTree = snpgdsCutTree(snpHCluster, z.threshold=15, outlier.n=5, n.perm = 5000, samp.group=NULL,
                         col.outlier="red", col.list=NULL, pch.outlier=4, pch.list=NULL,label.H=FALSE, 
                         label.Z=TRUE, verbose=TRUE)
snpgdsDrawTree(cutTree,type="z-score", main="Dataset 1")
snpgdsDrawTree(cutTree, main = "Dataset 1",edgePar=list(col=rgb(0.5,0.5,0.5,0.75),t.col="black"),
               y.label.kinship=T,leaflab="perpendicular",y.kinship.baseline = 0.8)

#pca
sample.id = read.gdsn(index.gdsn(genofile, "sample.id"))

pop_code = read.gdsn(index.gdsn(genofile, "sample.id"))
                      
pca = snpgdsPCA(genofile)
                      
tab = data.frame(sample.id = pca$sample.id,
                 pop = factor(pop_code)[match(pca$sample.id, sample.id)],
                 EV1 = pca$eigenvect[,1],
                 EV2 = pca$eigenvect[,2],stringsAsFactors = FALSE)

plot(tab$EV2, tab$EV1, col=as.integer(tab$pop),xlab="eigenvector 2", ylab="eigenvector 1")
legend("topleft", legend=levels(tab$pop), pch="o", col=1:nlevels(tab$pop))
