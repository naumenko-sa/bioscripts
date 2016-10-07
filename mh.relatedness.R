#calculate kinship for the families with SNPRelate

#example: https://www.biostars.org/p/83232/ - misleading
#tutorial: 
#http://www.bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html#format-conversion-from-vcf-files

#install SNPRelate as described here:
#http://www.bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html#installation-of-the-package-snprelate

library(gdsfmt)
library(SNPRelate)

setwd("~/Desktop/project_mh/2016-10-06_families_check/")

#prepare multisample vcf with bcftools merge 
#biallelic by default
snpgdsVCF2GDS("dataset1.vcf", "dataset1.gds")
snpgdsSummary("dataset1.gds")
genofile = snpgdsOpen("dataset1.gds")

#LD based SNP pruning
set.seed(1000)
snpset = snpgdsLDpruning(genofile,ld.threshold = 0.5)
snp.id=unlist(snpset)

dissMatrix  =  snpgdsIBS(genofile , sample.id=NULL, snp.id=snp.id, autosome.only=TRUE,remove.monosnp=TRUE, 
                         maf=NaN, missing.rate=NaN, num.thread=2, verbose=TRUE)
snpgdsClose(genofile)
snpHCluster =  snpgdsHCluster(dissMatrix, sample.id=NULL, need.mat=TRUE, hang=0.01)

cutTree = snpgdsCutTree(snpHCluster, z.threshold=15, outlier.n=5, n.perm = 5000, samp.group=NULL,
                         col.outlier="red", col.list=NULL, pch.outlier=4, pch.list=NULL,label.H=FALSE, 
                         label.Z=TRUE, verbose=TRUE)
snpgdsDrawTree(cutTree, main = "Dataset 1",edgePar=list(col=rgb(0.5,0.5,0.5,0.75),t.col="black"),
               y.label.kinship=T,leaflab="perpendicular")

