# calculates kinship for the families using SNPRelate

# example: https://www.biostars.org/p/83232/ - misleading
# tutorial: 
# http://www.bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html#format-conversion-from-vcf-files
# https://www.biostars.org/p/138694/
# docs:
# https://bioconductor.org/packages/release/bioc/manuals/SNPRelate/man/SNPRelate.pdf
# much simpler : vcftools --relatedness2
# https://www.biostars.org/p/111573/
# 1st degree ~0.25, and 2nd-degree ~0.125, and 3rd degree 0.0625.
# "Unrelated" parents can reach values as high as ~0.04 in my experience.

# install SNPRelate as described here:
# http://www.bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html#installation-of-the-package-snprelate

# Data preparation:
# 1.remove VEP annotations with vt
# 2.extract chrom 1-5 in the case of many samples (50)
# 3.extract SNPs only (problems with decomposed indels in bcftools merge)
# 3.prepare multisample vcf with bcftools merge  

installation <- function(){
    source("https://bioconductor.org/biocLite.R")
    biocLite("SNPRelate")
}

init <- function(){
    library(gdsfmt)
    library(SNPRelate)
}

tutorial <- function(){
    genofile <- snpgdsOpen(snpgdsExampleFileName())
    snpgdsSummary(snpgdsExampleFileName())
    pop_code <- read.gdsn(index.gdsn(genofile, path="sample.annot/pop.group"))
    table(pop_code)
    set.seed(1000)
    snpset <- snpgdsLDpruning(genofile, ld.threshold = 0.2)
    names(snpset)
    snpset.id <- unlist(snpset)
    
    pca <- snpgdsPCA(genofile, snp.id = snpset.id, num.thread = 2)
    pc.percent <- pca$varprop*100
    head(round(pc.percent,2))
    
    tab <- data.frame(sample.id = pca$sample.id,
                     EV1 = pca$eigenvect[,1],
                     EV2 = pca$eigenvect[,2],
                     stringsAsFactors = F)
    
    head(tab)
    plot(tab$EV2, tab$EV1, xlab = "eigenvector 2", ylab = "eigenvector 1")
    
    sample.id = read.gdsn(index.gdsn(genofile, "sample.id"))
    pop_code = read.gdsn(index.gdsn(genofile, "sample.annot/pop.group"))
    
    head(cbind(sample.id, pop_code))
    
    tab <- data.frame(sample.id = pca$sample.id,
                      pop = factor(pop_code)[match(pca$sample.id, sample.id)],
                      EV1 = pca$eigenvect[,1],    # the first eigenvector
                      EV2 = pca$eigenvect[,2],    # the second eigenvector
                      stringsAsFactors = F)
    head(tab)
    
    plot(tab$EV2, tab$EV1, col = as.integer(tab$pop), 
         xlab = "eigenvector 2", ylab = "eigenvector 1")
    legend("bottomright", legend = levels(tab$pop), pch = "o", col = 1:nlevels(tab$pop))
          
    #relatedness analysis
    sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
    YRI.id <- sample.id[pop_code == "YRI"]
    ibd <- snpgdsIBDMoM(genofile, sample.id = YRI.id, snp.id = snpset.id,
                        maf = 0.05, missing.rate = 0.05, num.thread = 2)
    
    ibd.coeff <- snpgdsIBDSelection(ibd)
    head(ibd.coeff)
    
    plot(ibd.coeff$k0, ibd.coeff$k1, xlim = c(0,1), ylim = c(0,1),
         xlab = "k0", ylab = "k1", main = "YRI samples (MoM)")
    lines(c(0,1), c(1,0), col = "red", lty = 2)
    
    
    dissMatrix <- snpgdsIBS(genofile, snp.id = snpset.id, autosome.only = T, 
                             maf = NaN, missing.rate = NaN, num.thread = 2, 
                             verbose=T, sample.id = YRI.id)
    
    snpHCluster <- snpgdsHCluster(dissMatrix, sample.id = NULL, need.mat = TRUE, hang = 0.01)
    
    #outlier.n=5
    cutTree <- snpgdsCutTree(snpHCluster, z.threshold = 15, outlier.n = 10, n.perm = 5000, 
                             samp.group = NULL, col.outlier = "red", col.list = NULL, 
                             pch.outlier = 4, pch.list = NULL, label.H = F, label.Z = T, 
                             verbose = T)
    snpgdsDrawTree(cutTree, type = "z-score")
    png(paste0(prefix, ".png"))
    snpgdsDrawTree(cutTree, main = "YRI", edgePar = list(col = "black", t.col = "black"),
                   y.label.kinship = T, leaflab = "perpendicular", outlier.n = 0)
    
    snpgdsClose(genofile)
}

# plots pca and dendrogram for samples listed in samples.txt in dataset.gds
# dataset.gds should be in the current directory: do coversion vcf -> gds first
plot_relatedness_picture <- function(samples.txt){
    #test:
    #samples.txt = "all.samples.txt"
    #samples.txt = "c4r_24.samples.txt"
    prefix <- gsub(".txt","",samples.txt,fixed = T)  
    samples <- unlist(read.table(samples.txt, stringsAsFactors = F))
    
    #biallelic by default
  
    #snpgdsSummary("dataset.gds")
    genofile <- snpgdsOpen("dataset2.gds")

    #LD based SNP pruning
    set.seed(1000)
    #default threshold is 0.2, for many samples it should be relaxed to 0.5
    snpset <- snpgdsLDpruning(genofile, ld.threshold = 0.2, sample.id = samples)
    snpset.id <- unlist(snpset)

    #pca
    pca <- snpgdsPCA(genofile, snp.id = snpset.id, num.thread = 2, sample.id = samples)

    pc.percent <- pca$varprop*100
    head(round(pc.percent,2),20)
    
    tab <- data.frame(sample.id = pca$sample.id,
                 EV1=pca$eigenvect[,1],
                 EV2=pca$eigenvect[,2],
                 stringsAsFactors = F
    )

    png(paste0(prefix, ".pca.png"))
    plot(tab$EV2,tab$EV1, xlab = "eigenvector 2", ylab = "eigenvector 1")
    text(tab$EV2,tab$EV1, tab$sample.id)
    dev.off()
    
    #family.id = c(1,1,1,2,2,2,3,3,3,4,5,6,7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,12,12,12,
    #              13,14,15,16,17,18,19,19,19,20,20,20,25,25,25,26,26,26,27,27,27,
    #              28,28,28,29,29,29,30,30,30)
    #family.id = c (1,1,2,2,2,3,3,4,5,6,6,6,7,7,7,8,9,9,9,10,10,11,11,12,13,13,14,14,14,15,16,17,17,17,18,18,19,19,19)
    family.id <- c(1,1,1)
                   
    ibd.robust <- snpgdsIBDKING(genofile, snp.id = snpset.id, num.thread = 2, 
                               family.id = family.id)
    dat <- snpgdsIBDSelection(ibd.robust, kinship.cutoff = 0.2)
    
    dissMatrix <- snpgdsIBS(genofile, snp.id = snpset.id, autosome.only = T,
                            remove.monosnp = F, maf = NaN, missing.rate = NaN, 
                            num.thread = 2, verbose = T, sample.id = samples)

    snpHCluster <- snpgdsHCluster(dissMatrix, sample.id = NULL, need.mat = TRUE, hang = 0.01)

    #outlier.n=5
    cutTree = snpgdsCutTree(snpHCluster, z.threshold = 5, outlier.n = 0, n.perm = 5000, 
                            samp.group = NULL, col.outlier = "red", col.list = NULL, 
                            pch.outlier = 4, pch.list = NULL, label.H = T, label.Z = T, verbose = T)
    
    snpgdsDrawTree(cutTree, type = "z-score")
    
    png(paste0(prefix,".png"), width = 1000)
    snpgdsDrawTree(cutTree, main = prefix, edgePar = list(col = "black", t.col = "black"),
               y.label.kinship = T, leaflab = "perpendicular", outlier.n = 0)
    dev.off()

    snpgdsClose(genofile)
}

init()
setwd("~/Desktop/work")
snpgdsVCF2GDS("C1A-106.vcf", "dataset2.gds")
plot_relatedness_picture("samples.txt")
