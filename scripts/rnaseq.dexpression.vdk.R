init = function()
{
  library(edgeR)
  setwd("~/Desktop/project_vdk/")
  source("~/Desktop/bioscripts/rnaseq.dexpression.R")
}

merge_counts()

samples = read.table("samples.txt", quote="\"", stringsAsFactors=F)
samples = samples[,1]

outliers = c("S19","S20","S21")

samples_no_outliers = setdiff(samples,outliers)
samples = samples_no_outliers


all_counts <- read.csv("~/Desktop/project_vdk/all_counts.txt", sep="", stringsAsFactors=F, row.names = 1)

n_samples = length(samples)
group=factor(c(rep(1,floor(n_samples/2)),rep(2,ceiling(n_samples/2))))
#patient = factor(c("511","511","511","523","523","523",
#                  "511","511","511","523","523","523"))

x=all_counts[samples]
y=DGEList(counts=x,group=group,genes=row.names(x),remove.zeros = T)

png("no_outliers.mds.png",width=2000)
plotMDS(y,cex=2)
dev.off()

#filter - 1 or 0.5
filter=1
keep=rowSums(cpm(y)>filter) >= n_samples/2
y=y[keep,,keep.lib.sizes=F]




y= calcNormFactors(y)

logcpm = cpm(x,prior.count=1,log=T)
