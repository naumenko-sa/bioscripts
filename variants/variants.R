# http://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html
installation <- function(){
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("VariantAnnotation", version = "3.8")
}

init <- function(){
    setwd("~/Desktop/work/uoft_tutorial_2019/")
    library(VariantAnnotation)
}

overview <- function(){
    file <- system.file("extdata", "chr22.vcf.gz", package = "VariantAnnotation")   
    vcf <- readVcf(file, "hg19")
    
    vcf 
    header <- header(vcf)
    
    samples(header(vcf))
    
    #genomic positions - genomic ranges
    length(rowRanges(vcf))
    gr <- rowRanges(vcf)
    # extract first 10 variants
    gr[1:10,c("QUAL","FILTER")]
    gr$QUAL
    df <- as.data.frame(values(gr))
    
    #access functions
    alt(vcf)
    qual(vcf)
    ref(vcf)
    
    #genotype data
    geno(vcf)
    genotypes <- as.data.frame(geno(vcf)$GT)
    nrow(genotypes)
    
    HG00096_HOM <- genotypes[genotypes$HG00096 == "1|1",]
    
    #info data = variant wise
    info(header(vcf))
    info_description <- info(header(vcf))
    
    subset = info(vcf)[1:4,1:5]
    subset["rs7410291",]
    info(vcf)["rs7410291","SVLEN"]
    
    region <- GRanges(seqnames="22", ranges = IRanges(start = c(50301422, 50989541),
                                                    end = c(50312106, 51001328),
                                                    names = c("gene1","gene2")))
    
    tab <- TabixFile(file)
    vcf.subset <- readVcf(tab, "hg19", param = region)
    genotypes <- geno(vcf.subset)$GT
    nrow(genotypes)
    
    writeVcf(vcf.subset,"subset.vcf")
    bgzip("subset.vcf", "subset.vcf.gz",overwrite=T)
    indexTabix("subset.vcf.gz",format="vcf")
    
    #don't be confused with write.vcf from bedr!
}

personal_genome <- function(){
    #vcf <- readVcf("PGPC_0001_S1.flt.vcf.gz","hg19")
 
    panel <- read.csv("~/cre/data/primary_immunodeficiency.csv")
    attach(panel)
    panel <- panel[,c("PanelAPP.EnsemblId.GRch37","PanelAPP.Gene.Symbol")]
    write.csv(panel, "immunodeficiency.csv", row.names = F)
    
    vcf_file <- "PGPC_0001_S1.flt.vcf.gz"
    
    panel_bed <- read.table("immunodeficiency.bed", stringsAsFactors = F)
    region <- GRanges(seqnames = panel_bed$V1, ranges = IRanges(start = panel_bed$V2,
                                                      end = panel_bed$V3,
                                                      names = panel_bed$V4))
    tab <- TabixFile(vcf_file)
    vcf_subset <- readVcf(tab, "hg19", param = region)
    
    samples(header(vcf_subset))
    
    writeVcf(vcf_subset,"PGPC_0001.subset.vcf")
    genotypes <- geno(vcf_subset)$GT
    nrow(genotypes)
    
    info(header(vcf_subset))
    
    #gr <- head(rowRanges(vcf_subset),100)
    gr <- rowRanges(vcf_subset)
    df <- data.frame(chr = seqnames(gr),
                     pos = start(gr))
    df <- cbind(df,gr$REF)
    df <- cbind(df,gr$ALT)
    df <- cbind(df,gr$QUAL)
    df <- cbind(df,gr$FILTER)
}