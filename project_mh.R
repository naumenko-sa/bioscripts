################################################
##  Malignant hyperthermia
################################################

mh_panel=c("CACNA1S","RYR1","STAC3","TRDN","ASPH","JPH2","CASQ1","ATP2A1","ATP2A2","CALM1","FKBP1A")

library(pheatmap)
library(RColorBrewer)
source("~/Desktop/bioscripts/rnaseq.load_rpkm_counts.R")

gene_lengths = read.delim("~/Desktop/project_muscular/reference/gene_lengths.txt", stringsAsFactors=F, row.names=1)
gtex_rpkm = read.csv("~/Desktop/project_muscular/reference/gtex.muscle_gene.rpkm", sep="", stringsAsFactors = F)

#also uses functions from dexpression.muscular.R

setwd("~/Desktop/project_mh/2017-02-01_RYR1_expression/")

samples = unlist(read.table("samples.txt", quote="\"", stringsAsFactors=F))

s1130_BD_B175 = load_rpkm_counts("1130-BD-B175.counts.rpkm")
s1258_AC_A79 = load_rpkm_counts("1258-AC-A79.counts.rpkm")
s1275_BK_B225 = load_rpkm_counts("1275-BK-B225.counts.rpkm")  
s1388_MJ_M219 = load_rpkm_counts("1388-MJ-M219.counts.rpkm")  
s2180_CB_C365 = load_rpkm_counts("2180-CB-C365.counts.rpkm")
s6087_PD_B317G = load_rpkm_counts("6087-PD-B317G.counts.rpkm")
s6100_DD_D360 = load_rpkm_counts("6100-DD-D360.counts.rpkm")


s1258_AC_A79$external_gene_name = NULL
s1275_BK_B225$external_gene_name = NULL
s1388_MJ_M219$external_gene_name = NULL
s2180_CB_C365$external_gene_name = NULL
s6087_PD_B317G$external_gene_name = NULL
s6100_DD_D360$external_gene_name = NULL

samples = merge_row_names(s1130_BD_B175,s1258_AC_A79)
samples = merge_row_names(samples,s1275_BK_B225)
samples = merge_row_names(samples,s1388_MJ_M219)
samples = merge_row_names(samples,s2180_CB_C365)
samples = merge_row_names(samples,s6087_PD_B317G)
samples = merge_row_names(samples,s6100_DD_D360)

write.table(samples,"mh.rpkms",quote=F,sep = "\t")
#sometimes contains duplicate entries, delete them
s62_AF_S5 = read.delim("62_AF_S5.rpkms.controls.txt", row.names=1, stringsAsFactors=F)

breaks = c(0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,2000)
plot_panel(mh_panel, samples, gtex_rpkm, "MH_genes_panel.png","Malignant hyperthermia genes RPKM",breaks)
