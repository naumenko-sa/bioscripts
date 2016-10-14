#calculate kinship for the families with SNPRelate

#example: https://www.biostars.org/p/83232/ - misleading
#tutorial: 
#http://www.bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html#format-conversion-from-vcf-files

#install SNPRelate as described here:
#http://www.bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html#installation-of-the-package-snprelate

library(gdsfmt)
library(SNPRelate)

#setwd("~/Desktop/project_mh/2016-10-06_families_check/")

setwd("~/Desktop/project_cheo/2016-10-14_family_check/")
#prepare multisample vcf with bcftools merge 
#biallelic by default
#snpgdsVCF2GDS("dataset2.vcf", "dataset2.gds")
snpgdsVCF2GDS("forge4.vcf", "dataset2.gds")
snpgdsSummary("dataset2.gds")
genofile = snpgdsOpen("dataset2.gds")

#LD based SNP pruning
set.seed(1000)
snpset = snpgdsLDpruning(genofile,ld.threshold = 0.5)
snp.id=unlist(snpset)

dissMatrix  =  snpgdsIBS(genofile , sample.id=NULL, snp.id=snp.id, autosome.only=TRUE,remove.monosnp=TRUE, 
                         maf=NaN, missing.rate=NaN, num.thread=2, verbose=TRUE)
snpgdsClose(genofile)
snpHCluster =  snpgdsHCluster(dissMatrix, sample.id=NULL, need.mat=TRUE, hang=0.01)

cutTree = snpgdsCutTree(snpHCluster, z.threshold=15, outlier.n=5, n.perm = 5000, samp.group=NULL,
                         col.outlier="red",col.list=NULL, pch.outlier=4, pch.list=NULL,label.H=FALSE, 
                         label.Z=TRUE, verbose=TRUE)
snpgdsDrawTree(cutTree, main = "Forge 4 samples",edgePar=list(col="black",t.col="black"),
               y.label.kinship=T,leaflab="perpendicular",outlier.n = 0)

