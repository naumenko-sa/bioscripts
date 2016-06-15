# https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
library(edgeR)
setwd("~/Desktop/project_katie_csc")
#1T4R2,1T4R3 - no 1T4T1
TPOINT=12
TR1=2
TR2=4

#SG511_ven_hi_2_26
#SG511_ven_hi_4_13
#SG511_ven_hi_4_27
#SG511_ven_lo_2_26
#SG511_ven_lo_4_13
#SG511_ven_lo_4_27
#SG523_ven_hi_2_27
#SG523_ven_hi_4_10
#SG523_ven_hi_4_24
#SG523_ven_lo_2_27
#SG523_ven_lo_4_10
#SG523_ven_lo_4_24


all_counts=read.delim("combined.counts",row.names="id")


#x=all_counts[c("SG511_ven_hi_2_26","SG511_ven_hi_4_13","SG511_ven_hi_4_27",
#               "SG511_ven_lo_2_26","SG511_ven_lo_4_13","SG511_ven_lo_4_27")]

#x=all_counts[c("SG523_ven_hi_2_27","SG523_ven_hi_4_10","SG523_ven_hi_4_24",
 #              "SG523_ven_lo_2_27","SG523_ven_lo_4_10","SG523_ven_lo_4_24")]

#x=all_counts[c("SG511_ven_hi_2_26","SG523_ven_hi_2_27",
#               "SG511_ven_lo_2_26","SG523_ven_lo_2_27")]

#x=all_counts[c("SG511_ven_hi_4_13","SG523_ven_hi_4_10",
#               "SG511_ven_lo_4_13","SG523_ven_lo_4_10")]

#x=all_counts[c("SG511_ven_hi_4_27","SG523_ven_hi_4_24",
#               "SG511_ven_lo_4_27","SG523_ven_lo_4_24")]

#x=all_counts[c("X1T3R1","X1T3R2","X1T3R3","X1T4R2","X1T4R3")]

x=all_counts[c("SG511_ven_hi_2_26","SG511_ven_hi_4_13","SG511_ven_hi_4_27",
"SG523_ven_hi_2_27","SG523_ven_hi_4_10","SG523_ven_hi_4_24",
               "SG511_ven_lo_2_26","SG511_ven_lo_4_13","SG511_ven_lo_4_27",
              "SG523_ven_lo_2_27","SG523_ven_lo_4_10","SG523_ven_lo_4_24")]


#group=factor(c(1,1,1,2,2,2))
#group=factor(c(1,1,2,2))
group=factor(c(1,1,1,1,1,1,2,2,2,2,2,2))

y=DGEList(counts=x,group=group)

y=calcNormFactors(y)

design=model.matrix(~group)

y=estimateDisp(y,design)

nfilename="all_hi_vs_lo.txt"
write.table(y$counts,nfilename)

fit=glmFit(y,design)
lrt=glmLRT(fit)
topTags(lrt,p.value=0.05,n=50)

efilename="all_hi_vs_lo.genes.txt"
write.table(topTags(lrt,p.value=0.05,n=50),efilename)
