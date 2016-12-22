###################################################################
############   Exon counts for Muscular project
#http://bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.pdf
#https://support.bioconductor.org/p/71345/
#JunctionSeq is better
#http://nar.oxfordjournals.org/content/early/2016/06/01/nar.gkw501.full
#http://stackoverflow.com/questions/21506724/how-to-plot-overlapping-ranges-with-ggplot2
#intron retention
#http://seqanswers.com/forums/showthread.php?t=42420
###################################################################
library(DEXSeq)
library(GenomicRanges)
library(IRanges)
library(ggplot2)
library(JunctionSeq)
library(edgeR)

init = function()
{
  library("biomaRt")  
  grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", 
                   path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
  datasets=listDatasets(grch37)
  
  grch37 = useDataset(grch37,dataset="hsapiens_gene_ensembl")
  
  attributes=listAttributes(grch37)
  filters=listFilters(grch37)
}
get_exons = function(ensembl_transcript_id)
{
  getBM(
    attributes=c('exon_chrom_start','exon_chrom_end','genomic_coding_start','genomic_coding_end'),
    filters=c('ensembl_transcript_id'),
    values=list(ensembl_transcript_id),
    mart=grch37)
}

setwd("~/Desktop/project_muscular/")

#remove headers from count files!
samples = c("fibroblast5","myotubes5","muscle5","muscle1","muscle2","1130-BD-B175",
            "1275-BK-B225","2180-CB-C365","6100-DD-D360",
            "1258-AC-A79","1388-MJ-M219","6087-PD-B317G")

#samples = c("muscle1","muscle2")
#condition = c("knockdown","control")

condition = rep("knockdown",12)
              #,rep("control",7))

sampleTable = data.frame(row.names=samples,
                         condition=condition)

files = paste0(samples,".dexseq")

t=DEXSeqDataSetFromHTSeq(countfiles=files,
                         sampleData = sampleTable, 
                         design =~ sample + exon,
                         flattenedfile = c("ref-transcripts.dexseq.gff3"))
#+condition:exon

#lama2, ryr1
genes = c("ENSG00000196569","ENSG00000196218","ENSG00000198947") 
gene_names = c("lama2","ryr1","dmd")

dxd = t[geneIDs(t) %in% genes[3],]
dxd <-estimateSizeFactors(dxd)

exonic_chunks_values = rowRanges(dxd)$exonBaseMean
exonic_chunks_ranges = ranges(rowRanges(dxd))
ir = data.frame(start=start(ranges(rowRanges(dxd))),end=end(ranges(rowRanges(dxd))))

exons = get_exons("ENST00000421865")

exons = get_exons("ENST00000355481")
exons = exons[,c("genomic_coding_start","genomic_coding_end")]
ir_exons = IRanges(exons$genomic_coding_start,exons$genomic_coding_end)
ol = findOverlaps(exonic_chunks_ranges,ir_exons,select="first")
overlap = as.matrix(ol)

png(paste0(gene_names[1],".barplot.png"),width=1000)
barplot(exonic_chunks_values,names.arg=overlap[,1],las=2,main=gene_names[1])

barplot(counts(dxd,normalized=T)[,c(1:5)],names.arg=overlap[,1],las=2,main=gene_names[1])

sample_n=5
plot(counts(dxd,normalized=T)[,4],
     type='p',main=gene_names[2],xlab="Exon",ylab="Expression",sub=samples[sample_n],ylim=c(0,200))
xaxt='n'
text(x=1:125, y=-200, overlap[,1], cex=0.8, srt=90, xpd=TRUE)

#125 in ryr1

cnts = as.data.frame(counts(dxd,normalized=T))
ggplot(as.data.frame(counts(dxd,normalized=T)[,c(1:2)]))

axis(1,at=1:73,labels=overlap[,1],cex.axis=0.7)


dev.off()

#colData(dxd)
rowRanges(dxd)[1]
#featureCounts(dxd)[10,]



#for additionalannotation
#gr1 = GRanges(c("6"),
#              IRanges(start=c(129204342,129371063),	
#                      end=c(129204502,129371233),names=c("1","2")))
#grl <- GRangesList("Exons" = gr1)

png(paste0(gene_names[3],".png"),width=3000)

plotDEXSeq(dxd, genes[3], legend=T, cex.axis=1.2, cex=1.3,
            lwd=2, displayTranscripts = T,
            expression = F, norCounts = T, names = T,
            color.samples=c("yellow","orange","red","green","blue",rep("gray",7)))


axis(3, at=exons$genomic_coding_start, labels=rownames(exons))


dev.off()

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


#JunctionSeq
#http://hartleys.github.io/JunctionSeq/doc/example-walkthrough.pdf
library(QoRTs)

#load data
path="/home/sergey/Desktop/project_muscular/jseq_muscle2"
setwd(path)
res = read.qc.results.data("",decoder.files="decoder.txt", calc.DESeq2 = F, calc.edgeR = T)
makeMultiPlot.all(res,outfile.dir = path,plot.device.name="pdf")
get.size.factors(res, outfile="sizeFactors.txt",sf.method = c("edgeR"))

#gene expression analysis with edgeR
suppressPackageStartupMessages(library(edgeR))
decoder.bySample = read.table("decoder.txt",header=T,stringsAsFactors=F)
#directory="countTables"
files = paste0(decoder.bySample$unique.ID,"/QC.geneCounts.formatted.for.DESeq.txt.gz")

countData <- lapply(files, function(f){
  ct <- read.table(f,header=F,stringsAsFactors=F)$V2;
  ct <- ct[1:(length(ct)-5)]
})

countMatrix <- do.call(cbind.data.frame,countData)
colnames(countMatrix) <- decoder.bySample$unique.ID
rownames(countMatrix) <- read.table(files[1],header=F,
                                    stringsAsFactors=F)$V1[1:nrow(countMatrix)]
group=factor(decoder.bySample$group.ID)

design <- model.matrix(~group)
y <- DGEList(counts = countMatrix, group = group);
y <- calcNormFactors(y)
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)

ensembl_w_description <- read.delim2("~/Desktop/reference_tables/ensembl_w_description.txt", stringsAsFactors=FALSE)
top_tags = as.data.frame(topTags(lrt,n=2000))
top_tags = merge(top_tags,ensembl_w_description,by.x = "row.names",by.y="ensembl_gene_id",all.x=T)
write.table(top_tags,"de_results.txt",quote=F,row.names=F,sep=";")

#differential splicing analysis
#first run : rnaseq.qorts.merge_novel_splices.sh
library(JunctionSeq)
decoder <- read.table("decoder.txt", header=T,stringsAsFactors=F)
countFiles <- paste0(decoder$sample.ID, "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")
#long step - 1h-3h
jscs <- runJunctionSeqAnalyses(
  sample.files = countFiles,
  sample.names = decoder$sample.ID,
  condition = decoder$group.ID,
  flat.gff.file = "withNovel.forJunctionSeq.gff.gz",
  nCores = 3,
  verbose=TRUE,
  debug.mode = TRUE
)

writeCompleteResults(jscs,outfile.prefix="results",save.jscs = T)

buildAllPlots(
  jscs=jscs,
  FDR.threshold = 0.01,
  outfile.prefix = "results",
  variance.plot = TRUE,
  ma.plot = TRUE,
  rawCounts.plot=TRUE,
  verbose = TRUE)

buildAllPlotsForGene(geneID = "ENSG00000196569",jscs=jscs)
