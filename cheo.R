#https://cran.r-project.org/web/packages/VennDiagram/VennDiagram.pdf
setwd("~/Desktop/project_cheo/variant_calls/omim")
setwd("~/Desktop/project_cheo/variant_calls/omim/snps")
setwd("~/Desktop/project_cheo/variant_calls/omim/indels")
library("VennDiagram")
venn.plot=draw.pairwise.venn(945+4958,132+4958,4958,
              c("BCBIO","Jacek"),fill=c("blue","red"),lty="blank",
              cex = 2, cat.cex=2, cat.just = list(c(-1,-1),c(1,1)),
              ext.length = 0.3,  ext.line.lwd=2,
              ext.text = F, main="Sdf"
              )
venn.plot=draw.pairwise.venn(46807953, 54328695, 44991414,
                             c("Nimblegen.capture","Agilent"),
                             fill=c("blue","red"),cex=2,cat.cex=2)

venn.plot3=draw.triple.venn(54328695,46807953,65824553,
                            39517373,43777891,44991414,
                            37679751,
                            c("Agilent","Nimblegen.capture","Nimblegen.empirical")
)

venn.plot3=draw.triple.venn(1,2,3,12,23,13,123,
  c("Agilent","Nimblegen.capture","Nimblegen.empirical")
)

#example
for (pipeline in c("genap","bcbio","jacek","cpipe"))
{
  for (sample in c("S1","S2","S4","S5","S6","S7","S8"))
  {
    agilent=unlist(read.table(paste0("S1.",pipeline,".agilent.omim.variants.indels")))
    nimblegen=unlist(read.table(paste0("S1.",pipeline,".nimblegen.omim.variants.indels")))
    #illumina=unlist(read.table(paste0("S1.",pipeline,".illumina.txt")))

    x=list(agilent,nimblegen)
    names=list("Agilent","Nimblegen")
    png(paste0(sample,".",pipeline,".png"))
    grid.draw(venn.diagram(x,NULL,category.names=names))
    dev.off()
  }
}

for (platform in c("agilent","nimblegen"))
{
  for (sample in c("S1","S2","S4","S5","S6","S7","S8"))
  {
    
  #sample="S3"
  #platform="illumina"
  jacek=unlist(read.table(paste0(sample,".jacek.",platform,".omim.variants.indels")))
  cpipe=unlist(read.table(paste0(sample,".cpipe.",platform,".omim.variants.indels")))
  bcbio=unlist(read.table(paste0(sample,".bcbio.",platform,".omim.variants.indels")))
  genap=unlist(read.table(paste0(sample,".genap.",platform,".omim.variants.indels")))

  x=list(jacek,cpipe,bcbio,genap)
  names=list("jacek","cpipe","bcbio","genap")

  png(paste0(sample,".",platform,".png"))
  grid.draw(venn.diagram(x,NULL,category.names=names))
  dev.off()
  }
}

# intersect bed files
# http://davetang.org/muse/2013/01/02/iranges-and-genomicranges/
source("http://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
library("GenomicRanges")
setwd("coverage/nimblegen")
capture=read.table("nimblegen.capture",header=F)
colnames(capture)=c('chr','start','end')
capture.bed = with(capture, GRanges(chr, IRanges(start+1, end)))

empirical=read.table("nimblegen.empirical",header=F)
colnames(empirical)=c('chr','start','end')
empirical.bed = with(capture, GRanges(chr, IRanges(start+1, end)))

omim=read.table("omim.orphanet.goodnames.v2.bed",header=F)
colnames(empirical)=c('chr','start','end')
omim.bed = with(capture, GRanges(chr, IRanges(start+1, end)))

bed.intersect = intersect(omim.bed,capture.bed)


#bed files
setwd("venn_diagrams/bed_intersection/")
omim=unlist(read.table("omim.orphanet.goodnames.v2.bed.exons"))
agilent=unlist(read.table("omim_vs_agilent.50percent.wo.exons"))
nimblegen=unlist(read.table("omim_vs_nimblegen.capture.50percent.wo.exons"))
illumina=unlist(read.table("omim_vs_illumina.50percent.wo.exons"))

x=list(omim,agilent,nimblegen,illumina)
names=list("omim","agilent","nimblegen","illumina")

grid.draw(venn.diagram(x,NULL,category.names=names))

#variants parameters
type="snps"
type="indels"
for (platform in c("agilent","nimblegen"))
{
  for (sample in c("S1","S2","S4","S5","S6","S7","S8"))
  {
    
    sample="S1"
    platform="agilent"
    cpipe.hom=read.table(paste0(sample,".cpipe.",platform,".omim.",type,".hom.AD"))
    cpipe.het=read.table(paste0(sample,".cpipe.",platform,".omim.",type,".het.AD"))
    
    bcbio.hom=read.table(paste0(sample,".bcbio.",platform,".omim.",type,".hom.AD"))
    bcbio.het=read.table(paste0(sample,".bcbio.",platform,".omim.",type,".het.AD"))
    
    genap.hom=read.table(paste0(sample,".genap.",platform,".omim.",type,".hom.AD"))
    genap.het=read.table(paste0(sample,".genap.",platform,".omim.",type,".het.AD"))

    genap.only = read.table("S1.genap.only.recode.vcf.AD")
    v=c(genap.hom,genap.het,genap.only)
    names=c("genap.hom","genap.het","genap.only")
    
    #indels
    v=c(cpipe.het,bcbio.het,genap.het)
    names=c("cpipe.het","bcbio.het","genap.het")
    
    png(paste0(sample,".",platform,".",type,".png"),width=1000)
    boxplot(v,names=names)
    dev.off()
  }
}

setwd("~/cluster/dorin_test")
b1=read.table("UNIQUE_to_BCBIO_1496461.vcf.AD")
b2=read.table("UNIQUE_to_BCBIO_1496462.vcf.AD")
b3=read.table("UNIQUE_to_BCBIO_1496463.vcf.AD")
g1=read.table("UNIQUE_to_DNASEQ_1496461.vcf.AD")
g2=read.table("UNIQUE_to_DNASEQ_1496462.vcf.AD")
g3=read.table("UNIQUE_to_DNASEQ_1496463.vcf.AD")

v=c(b1,b2,b3,g1,g2,g3)
names=c("b1","b2","b3","g1","g2","g3")

png(paste0(sample,".",platform,".",type,".png"),width=1000)
boxplot(v,names=names)

get_muscle_genes_coordinates = function()
{
    setwd("~/Desktop/project_cheo/")
    gene_list = unlist(read.table("cheo.muscular_genes", stringsAsFactors=F))
    setwd("~/Desktop/reference_tables/")
    all_exons = read.delim("protein_coding_genes.exons.bed", header=F, stringsAsFactors=F)
    muscular_exons = all_exons[all_exons$V4 %in% gene_list,]
    write.table(muscular_exons,"cheo.muscular_exons.bed",sep="\t",quote=F,row.names=F,col.names=F)
}

coverage.muscle_genes = function ()
{
    setwd("~/Desktop/project_cheo/2017-04-05_NextSeq_coverage/")
    samples = unlist(read.table("samples.txt", stringsAsFactors=F))
    
    coverage = read.delim(paste0(samples[1],".coverage"),header=T,stringsAsFactors = F)
    coverage = coverage[,c("gene","mean")]
    colnames(coverage)[2]=samples[1]
    
    
    for (sample in tail(samples,-1))
    {
        sample_coverage = read.delim(paste0(sample,".coverage"),header=T,stringsAsFactors = F)
        sample_coverage = sample_coverage[,c("gene","mean")]
        colnames(sample_coverage)[2]=sample
        coverage = cbind(coverage,sample_coverage[2])
    }
    row.names(coverage) = coverage$gene
    coverage$gene=NULL
    
    png("cheo.muscular_genes.coverage.part1.png",res=300,width=5000,height=2000)
    boxplot(t(coverage[1:60,]),las=2,cex.axis=0.8,main="Coverage in 60 samples of NextSeq for muscle gene panel,part 1")
    dev.off()
    
    png("cheo.muscular_genes.coverage.part2.png",res=300,width=5000,height=2000)
    boxplot(t(coverage[61:118,]),las=2,cex.axis=0.8,main="Coverage in 60 samples of NextSeq for muscle gene panel,part 2")
    dev.off()
}

#when looking at all genes, some samples may have no coverage
# CA0229.coverage 15325
# CA0246.coverage 17268
# CH0317.coverage 15833
# GM15262.coverage 14435
coverage.all_genes = function ()
{
    library("matrixStats")
    setwd("~/Desktop/project_cheo/2017-04-05_NextSeq_coverage/all_genes/")
    samples = unlist(read.table("samples.txt", stringsAsFactors=F))
  
    #hopefully 1st sample has most genes
    coverage = read.delim(paste0(samples[1],".coverage"),header=T,stringsAsFactors = F)
    coverage = coverage[,c("gene","mean")]
    colnames(coverage)[2]=samples[1]
    row.names(coverage)=coverage$gene
    coverage$gene=NULL
  
    for (sample in tail(samples,-1))
    {
        sample_coverage = read.delim(paste0(sample,".coverage"),header=T,stringsAsFactors = F)
        sample_coverage = sample_coverage[,c("gene","mean")]
        colnames(sample_coverage)[2]=sample
        coverage = merge(coverage,sample_coverage,by.x="row.names",by.y="gene",all.x=T)
        row.names(coverage) = coverage$Row.names
        coverage$Row.names=NULL
    }
  
    coverage[is.na(coverage)]=0  
    coverage$Mean = rowMeans(coverage)
    png("cheo.all_genes.coverage.png",res=300,width=5000,height=2000)
    boxplot(coverage,las=2,cex.axis=0.8,main="Coverage in 60 samples of NextSeq across all protein coding genes, no outliers",outline = F)
    dev.off()
    
    
    meds=rbind(colnames(coverage),colMedians(as.matrix(coverage)))
    write.table(meds,"medians.txt",col.names = F, quote = F,row.names = F)
}

omim_table_manipulation = function()
{
    setwd("~/Desktop/reference_tables/OMIM_2017-04-13/")
    mimTitles.percent <- read.delim2("mimTitles.percent.txt", comment.char="#")    
    genemap2 <- read.delim("genemap2.txt", comment.char="#")
    mimTitles.percent = merge(mimTitles.percent,genemap2, by.x="Mim.Number",by.y="Mim.Number",all.x=T,all.y=F)
    write.csv(mimTitles.percent,"mimTitles.percent.location.csv",row.names = F)
}

read_length_distribution = function(family)
{
    
    family_data = read.delim(paste0(family,".tsv"), stringsAsFactors=F)
    samples = unique(family_data$sample)
    
    for (sample in samples)
    {
        tmp = subset(family_data,sample==sample)
        tmp$sample=NULL
        print(paste0(sample," ",round(sum(tmp$Length*tmp$Count) / sum(tmp$Count))),quote=F)
    }
}

read_lengths = function()
{
    setwd("~/Desktop/project_cheo/2017-04-12_read_lengths/")
    families = unlist(read.table("projects.txt", stringsAsFactors=F))
    for (family in families)
    {
      read_length_distribution(family)
    }

    read_lengths <- read.csv("read_lengths.txt", sep="", stringsAsFactors=FALSE)
    read_lengths$id=NULL

    png("read_lengths_all_samples.png",width=2000)
    barplot(read_lengths$average_read_length, names.arg = read_lengths$sample,main = "Average read lengths for NextSeq samples is 134",
        las=2)
    dev.off()
}
