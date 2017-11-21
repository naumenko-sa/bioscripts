installation = function()
{
    source("https://bioconductor.org/biocLite.R")
    biocLite("VariantAnnotation")
}

init = function()
{
    setwd("~/Desktop/work")
    library(VariantAnnotation)
}

overview = function()
{
    file = system.file("extdata","chr22.vcf.gz",package="VariantAnnotation")   
    vcf = readVcf(file,"hg19")
    
    vcf 
    header=header(vcf)
    
    samples(header(vcf))
    
    #genomic positions
    length(rowRanges(vcf))
    
    alt(vcf)
    qual(vcf)
    ref(vcf)
    
    #genotype data
    geno(vcf)
    genotypes = geno(vcf)$GT
    
    #info data = variant wise
    info(header(vcf))
    info_description = info(header(vcf))
    
    subset = info(vcf)[1:4,1:5]
    subset["rs7410291",]
    info(vcf)["rs7410291","SVLEN"]
    
    region = GRanges(seqnames="22",ranges = IRanges(start=c(50301422,50989541),
                                                    end=c(50312106,51001328),
                                                    names = c("gene1","gene2")))
    
    tab=TabixFile(file)
    vcf.subset = readVcf(tab,"hg19",param=region)
    rowRanges(vcf.subset)
    genotypes = geno(vcf.subset)$GT
    
    writeVcf(vcf.subset,"subset.vcf")
    bgzip("subset.vcf", "subset.vcf.gz",overwrite=T)
    indexTabix("subset.vcf.gz",format="vcf")
    
    #don't be confused with write.vcf from bedr!
}

ashkenazim = function ()
{
    vcf = readVcf("Ashkenazim.vcf.gz","Grch37")
    samples(header(vcf))
}