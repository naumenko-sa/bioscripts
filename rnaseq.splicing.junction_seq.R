#JunctionSeq
#http://hartleys.github.io/JunctionSeq/doc/example-walkthrough.pdf
library(QoRTs)

#load data
#path="~/work/project_muscular/jseq_fibro5"
#setwd(path)
path=getwd()
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
#long step - 1h
jscs <- runJunctionSeqAnalyses(
  sample.files = countFiles,
  sample.names = decoder$sample.ID,
  condition = decoder$group.ID,
  flat.gff.file = "withNovel.forJunctionSeq.gff.gz",
  nCores = 1,
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
  rawCounts.plot=F,
  without.TX = F,
  expr.plot = F,
  normCounts.plot = T,
  plot.novel.junction.results = T,
  verbose = TRUE)


