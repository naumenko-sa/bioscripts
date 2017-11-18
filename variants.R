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