#calculates kinship for the families using SNPRelate

#example: https://www.biostars.org/p/83232/ - misleading
#tutorial: 
#http://www.bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html#format-conversion-from-vcf-files
#https://www.biostars.org/p/138694/

#install SNPRelate as described here:
#http://www.bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html#installation-of-the-package-snprelate


# 1.remove VEP annotations with vt
# 2.extract chrom 1-5
# 3.prepare multisample vcf with bcftools merge  

library(gdsfmt)
library(SNPRelate)

setwd("~/Desktop/project_cheo/2017-04-06_NextSeq_kinship/")

#biallelic by default
#snpgdsVCF2GDS("dataset2.vcf", "dataset2.gds")
#snpgdsVCF2GDS("S700.vcf", "dataset.gds")
snpgdsVCF2GDS("CHEO_0001.no_vep.vcf.gz", "dataset1.gds")
snpgdsSummary("dataset1.gds")
genofile = snpgdsOpen("dataset1.gds")

#LD based SNP pruning
set.seed(1000)
#default threshold is 0.2, for many samples it should be relaxed to 0.5
snpset = snpgdsLDpruning(genofile,ld.threshold = 0.2)
snp.id=unlist(snpset)

#pca
pca = snpgdsPCA(genofile,snp.id=snp.id,num.thread=2)

tab = data.frame(sample.id = pca$sample.id,
                 EV1=pca$eigenvect[,1],
                 EV2=pca$eigenvect[,2],
                 stringsAsFactors = F
)
plot(tab$EV2,tab$EV1,xlab="eigenvector 2",ylab="eigenvector 1")

#problematic: 
#sample.id=c("1419-MJ-M09","1391-MD-M09","1330-DH-C155","1284-MF-M116")
sample.id=c("57_SF","59_FF","60_PF","62_AF")
#outliers: 
#  "1337-UP-U03",
#  "2118-BS-B440" = 1418-KS-B119
#   1365-SD-S" = 1393-CS-D255
# "

sample.id=c("1333-SL-S700","1422-SE-S700",
            "1418-KS-B199","2115-BK-B199",
            "2215-DY-G358","2336-HC-G358",
            "1309-JZ-Z03",
            "1359-DW-B440",
            "1311-DV-D133",
            "1325-WM-W112",
            "1334-HJ-H203",
            "1336-HJ-H16",
            "1358-MT-M128",
            "1362-MD-M13",
            "1387-BA-B224",
            "1393-CS-D255",
            "1414-FJ-F129")


dissMatrix  =  snpgdsIBS(genofile, snp.id=NULL, autosome.only=TRUE,remove.monosnp=TRUE, 
                         maf=NaN, missing.rate=NaN, num.thread=2, verbose=TRUE)

#snpgdsClose(genofile)
snpHCluster =  snpgdsHCluster(dissMatrix, sample.id=NULL, need.mat=TRUE, hang=0.01)

cutTree = snpgdsCutTree(snpHCluster, z.threshold=15, outlier.n=5, n.perm = 5000, samp.group=NULL,
                         col.outlier="red",col.list=NULL, pch.outlier=4, pch.list=NULL,label.H=FALSE, 
                         label.Z=TRUE, verbose=TRUE)
snpgdsDrawTree(cutTree, main = "Varese, chr1_5",edgePar=list(col="black",t.col="black"),
               y.label.kinship=T,leaflab="perpendicular",outlier.n = 0)

snpgdsClose(genofile)
