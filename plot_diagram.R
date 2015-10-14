#run: Rscript plot_diagramR file.dat
args=commandArgs(trailingOnly=TRUE)
mydata=read.table(args[1])
library(lattice);
png(paste(args[1],".png",sep=""));
limits=c(min(mydata),max(mydata));
histogram(mydata$V1,xlab="Insert size",main=paste("Distribution of insert sizes,lib=",args[1],sep=""),endpoints=limits,xlim=limits);
dev.off();