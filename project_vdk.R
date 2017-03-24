setwd("~/Desktop/project_vdk/")

get_data = function()
{

    samples = read.table("samples.txt", quote="\"", stringsAsFactors=F)
    samples = samples[,1]
  
    for (sample in samples)
    {
        table_name = sample
        file_name = paste0(sample,".counts")
        assign(table_name,read.delim(file_name, stringsAsFactors=F))
    }
  
    samples.data = get(head(samples[1]))
  
    for (sample in tail(samples,-1))
    {
        table_name = sample
        table_itself = get(table_name)
        table_itself$symbol=NULL
        samples.data = merge(samples.data,get(table_name))
    }

    write.table(samples.data,"all_counts.txt",quote=F,row.names = F)
    
}

init = function()
{
  library(edgeR)
}

samples20 = c("S10","S11","S1","S12","S13","S14","S15","S16","S17",
            "S19","S20","S21","S2","S3","S4","S5","S6","S7","S8","S9")

samples = samples20

samples_no_outliers = c("S10","S11","S1","S12","S13","S14","S15","S16","S17",
              "S2","S3","S4","S5","S6","S7","S8","S9")

samples = samples_no_outliers

all_counts <- read.csv("~/Desktop/project_vdk/all_counts.txt", sep="", stringsAsFactors=F, row.names = 1)

n_samples = length(samples)
group=factor(c(rep(1,floor(n_samples/2)),rep(2,ceiling(n_samples/2))))
#patient = factor(c("511","511","511","523","523","523",
#                  "511","511","511","523","523","523"))

x=all_counts[samples]
y=DGEList(counts=x,group=group,genes=row.names(x),remove.zeros = T)

plotMDS(y,cex=.2)
#filter - 1 or 0.5
filter=1
keep=rowSums(cpm(y)>filter) >= n_samples/2
y=y[keep,,keep.lib.sizes=F]


png("no_outliers.mds.png",width=2000)
plotMDS(y,cex=2)
dev.off()

y= calcNormFactors(y)

logcpm = cpm(x,prior.count=1,log=T)
