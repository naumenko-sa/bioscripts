################################################
##  Malignant hyperthermia
################################################

library(pheatmap)
library(RColorBrewer)
source("~/Desktop/bioscripts/rnaseq.load_rpkm_counts.R")
source("~/Desktop/bioscripts/rnaseq.muscular_gene_panels.R")

setwd("~/Desktop/project_mh/2017-02-01_RYR1_expression/")

samples = unlist(read.table("samples.txt", quote="\"", stringsAsFactors=F))

plot_rpkm_table = function()
{
    s1130_BD_B175 = load_rpkm_counts("1130-BD-B175.counts.rpkm")
    s1258_AC_A79 = load_rpkm_counts("1258-AC-A79.counts.rpkm")  
    s1275_BK_B225 = load_rpkm_counts("1275-BK-B225.counts.rpkm")  
    s1388_MJ_M219 = load_rpkm_counts("1388-MJ-M219.counts.rpkm")  
    s2180_CB_C365 = load_rpkm_counts("2180-CB-C365.counts.rpkm")
    s6087_PD_B317G = load_rpkm_counts("6087-PD-B317G.counts.rpkm")
    s6100_DD_D360 = load_rpkm_counts("6100-DD-D360.counts.rpkm")


    for(sample in c(s1258_AC_A79,s1275_BK_B225,s1388_MJ_M219,s2180_CB_C365,s6087_PD_B317G,s6100_DD_D360))
    {
        sample$external_gene_name = NULL
    }

    counts = merge_row_names(s1130_BD_B175,s1258_AC_A79)
    for (next_sample in c(s1275_BK_B225,s1388_MJ_M219,s2180_CB_C365, s6087_PD_B317G, s6100_DD_D360))
    {
        counts = merge_row_names(counts,next_sample)
    }

    write.table(counts,"mh.rpkms",quote=F,sep = "\t")
    #sometimes contains duplicate entries, delete them
}

counts = read.delim("mh.rpkms", row.names=1, stringsAsFactors=F)

breaks = c(0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,2000)
plot_panel(mh_panel, samples, gtex_rpkm, "MH_genes_panel.png","Malignant hyperthermia genes RPKM",breaks)

plot_all_panels(counts,gtex_rpkm)
