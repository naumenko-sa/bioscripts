###################################################################
############   Exon counts for Muscular project
#http://bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.pdf
#https://support.bioconductor.org/p/71345/
###################################################################
library(DEXSeq)

setwd("~/Desktop/project_muscular/")

#remove headers from count files!
samples = c("fibroblast5","myotubes5","muscle5","muscle1","muscle2",
            "1130-BD-B175","1275-BK-B225","2180-CB-C365","6100-DD-D360",
            "1258-AC-A79","1388-MJ-M219","6087-PD-B317G")

samples = c("muscle1","muscle2")
condition = c("knockdown","control")

condition = c(rep("knockdown",5),rep("control",7))

sampleTable = data.frame(row.names=samples,
                         condition=condition)

files = paste0(samples,".dexseq")

t=DEXSeqDataSetFromHTSeq(countfiles=files,
                         sampleData = sampleTable, 
                         design =~ sample + exon + condition:exon,
                         flattenedfile = c("ref-transcripts.dexseq.gff3"))

#lama2, ryr1
genes = c("ENSG00000196569","ENSG00000196218") 

dxd = t[geneIDs(t) %in% genes,]

colData(dxd)

rowRanges(dxd)[45,]
featureCounts(dxd)[10,]

dxd = estimateSizeFactors(dxd)

plotDEXSeq( dxr1, genes[1], legend=TRUE, cex.axis=1.2, cex=1.3,
            lwd=2, displayTranscripts = T,
            expression = F, norCounts = T, names = T)
            color.samples=c("yellow","orange","red","green","blue",rep("gray",7)))


dxd = estimateDispersions(dxd)
plotDispEsts(dxd)
dxd = testForDEU(dxd)
dxd = estimateExonFoldChanges(dxd,fitExpToVar = "condition")
dxr1 = DEXSeqResults( dxd )

plotDEXSeq( dxr1, genes[1], legend=TRUE, cex.axis=1.2, cex=1.3,
            lwd=2, displayTranscripts = T,
            expression = F, norCounts = T, names = T,
            color.samples=c("yellow","orange","red","green","blue",rep("gray",7)))

formulaFullModel    =  ~ sample + exon + condition:exon
formulaReducedModel =  ~ sample + exon

dxd = estimateDispersions( dxd, formula = formulaFullModel )

dxd = testForDEU( dxd,
                  reducedModel = formulaReducedModel,
                  fullModel = formulaFullModel )

dxr2 = DEXSeqResults( dxd )

plotDEXSeq( dxr2, genes[1], legend=TRUE, cex.axis=1.2, cex=1.3,
            lwd=2, displayTranscripts = T,
            expression = T, norCounts = F, names = T)

