#all_counts = read.delim("annotated_combined.counts", row.names=1, stringsAsFactors=FALSE)
all_counts = read.csv("~/Desktop/project_muscular/all_counts.txt", sep="", stringsAsFactors=FALSE)

symbol = all_counts[,c(1,8)]
#all_counts = all_counts[,c(1,2,3,4,5,6,7,9,10,11,12,13,14,15)]

#path = "/counts/mh_filtered/"
get_data = function (samples_file, path)
{
    samples = read.table(samples_file, quote="\"", stringsAsFactors=F)

    for (sample in samples[,1])
    {
        table_name = paste0(sample,".filtered")
        file_name = paste0(getwd(),path,sample,".filtered.counts")
        assign(table_name,read.delim(file_name, stringsAsFactors=F))
    }

    samples.data = get(paste0(head(samples,1),".filtered"))

    for (sample in tail(samples,-1)[,1])
    {
        table_name = paste0(sample,".filtered")
        table_itself = get(table_name)
        table_itself$symbol=NULL
        samples.data = merge(samples.data,get(table_name))
    }
    return(samples.data)
}

mh.samples = get_data("mh.samples.txt","/counts/mh_filtered/")
muscular.samples = get_data("muscular.samples.txt","/counts/muscular_filtered/")
all.samples.data = merge(mh.samples,muscular.samples)

all.samples.data$id=NULL
all.samples.data$symbol=NULL
all.samples.total_counts = colSums(all.samples.data)
all.samples.non_zero = colSums(all.samples.data!=0)

png("Fig2.total_counts_after_filtration.png",width=1000)
op <- par(mar = c(10,4,4,2) + 0.1)
barplot(all.samples.total_counts,las=2,main="Total counts after filtration")
par(op)
dev.off()

png("Fig3.Genes_with_reads_after_filtration.png",width=1000)
op <- par(mar = c(10,4,4,2) + 0.1)
barplot(all.samples.non_zero,las=2,main="Genes with reads after filtration")
par(op)
dev.off()

bars = list()
for (col_name in colnames(all.samples.data))
{
    #col_name = 'X1258.AC.A79'
    col_data = all.samples.data[,col_name]
    col_data = log2(col_data[col_data!=0])
    bars = c(bars, list(col_data))
    #boxplot(col_data,las=2)
    #assign(paste0("bars$",col_name),col_data)
}

png("Fig4.Log2_counts_for_covered_genes.png",width=1000)
op <- par(mar = c(10,4,4,2) + 0.1)
boxplot(bars,names = colnames(all.samples.data),las=2,main="Log2_counts_for_covered_genes")
par(op)
dev.off()


norm_counts = cbind(symbol,cpm(all_counts[,c(2,3,4,5,6,7,9,10,11,12,13,14,15)]))
write.table(norm_counts,"norm_counts_all.txt",col.names=NA,quote=F)  

gene_locus = read.delim("genes_locus.txt",stringsAsFactors = F,header=T)
gene_panel = read.delim("genes_muscular.txt",stringsAsFactors = F,header=T)

muscular.locus = muscular.samples[muscular.samples$symbol %in% unlist(gene_locus),]
mh.locus = mh.samples[mh.samples$symbol %in% unlist(gene_locus),]
muscular.locus$id=NULL
muscular.locus$symbol=NULL
mh.locus$id=NULL
mh.locus$symbol=NULL

png("Fig5.Log2_counts_for_locus_mh_samples.png",width=1000)
boxplot(log2(mh.locus+1),main="Log2 counts for locus mh samples")
dev.off()

png("Fig6.Log2_counts_for_locus_muscular_samples.png",width=1000)
boxplot(log2(muscular.locus+1))
dev.off()

muscular.panel = muscular.samples[muscular.samples$symbol %in% unlist(gene_panel),]
mh.panel = mh.samples[mh.samples$symbol %in% unlist(gene_panel),]
muscular.panel$id=NULL
muscular.panel$symbol=NULL
mh.panel$id=NULL
mh.panel$symbol=NULL

png("Fig7.Log2_counts_for_panel_mh_samples.png",width=1000)
boxplot(log2(mh.panel+1),main="Log2 counts for locus mh samples")
dev.off()

png("Fig8.Log2_counts_for_panel_muscular_samples.png",width=1000)
boxplot(log2(muscular.panel+1))
dev.off()

work_counts=muscular_genes
work_counts = all_counts

#exploration
row.names(work_counts) = work_counts$symbol
work_counts$id = NULL
work_counts$symbol = NULL
total_counts = colSums(work_counts)
barplot(total_counts)
log_counts = log2(work_counts)
png("all_genes_counts.png")
par(mar=c(10,3,1,1))
boxplot(log_counts,las=2)
dev.off()

attach(work_counts)

x=muscular_genes[c("muscle1","X1130.BD.B175")]

x=muscular_genes[c(2,3,4,5,6,7,9,10,11,12,13,14,15)]
group=factor(c(1,2,3,4,5,6,7,8,9,10,11,12,13))

group=factor(c(1,2))

y=DGEList(counts=x,group=group)


y=calcNormFactors(y)

plotMDS(y)
bcv=0.2
#from A.thaliana experiment
dispersion=0.04
  
  design=model.matrix(~group)
  fit=glmFit(y,design,dispersion)
  lrt=glmLRT(fit,coef=2)

write.table(merge(work_counts,lrt$table,by="row.names"),"muscle1_vs_1130_panel.txt",col.names=NA,quote=F)  

expression_unfiltered = function()
{
    setwd("~/Desktop/project_muscular/counts/muscular_unfiltered/")
    annotated_combined <- read.delim("~/Desktop/project_muscular/counts/muscular_unfiltered/annotated_combined.counts", row.names=1, stringsAsFactors=FALSE)
    all_genes_lengths <- read.delim("~/Desktop/project_muscular/counts/muscular_unfiltered/all_genes_lengths", row.names=1, stringsAsFactors=FALSE)
    
    samples = c("myotubes5","fibroblast5")  
    counts = annotated_combined[,samples]

    counts = merge(counts,all_genes_lengths,by.x="row.names",by.y="row.names")   
    row.names(counts)=counts$Row.names
    counts$Row.names = NULL
    
    x=counts[samples]
    group = factor(c(1,1))
    y=DGEList(counts=x,group=group,remove.zeros = F)
    plotMDS(y)
    #generate counts in the other way
    rpkms = rpkm(y,counts$Length)
    
    rpkms = merge(rpkms,ensembl_w_description,by.x="row.names",by.y="row.names")
    row.names(rpkms)=rpkms$Row.names
    rpkms$Row.names=NULL
    
    panel.rpkm = rpkms[rpkms$external_gene_name %in% congenital_myopathy,]
    row.names(panel.rpkm) = panel.rpkm$external_gene_name
    panel.rpkm$external_gene_name=NULL
    panel.rpkm$Gene_description=NULL
    
    png("congenital_myopaties.png",res=100,width=1000)
    pheatmap(panel.rpkm,treeheight_row=0,treeheight_col=0,cellwidth = 40,
             display_number =T,cluster_rows=T, cluster_cols=T,
             main="Congenital myopaties, RPKM/2, unfiltered",breaks = c(0,10,50,100,500,1000,5000,10000,15000),
             colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(8))
    dev.off()
}

# plots the heatmap of title.png for gene_panel using sample_rpkm and gtex_rpkm
plot_panel= function(gene_panel, sample_rpkm, gtex_rpkm, filename,title, breaks)
{
    panel_rpkm = sample_rpkm[sample_rpkm$external_gene_name %in% gene_panel,]
    gtex_panel_rpkm = gtex_rpkm[gtex_rpkm$gene_name %in% gene_panel,]
  
    all_rpkm = merge(panel_rpkm,gtex_panel_rpkm,by.x='external_gene_name',by.y='gene_name')
    row.names(all_rpkm) = all_rpkm$external_gene_name
    all_rpkm$external_gene_name=NULL
  
    png(filename,res=100,width=600)
    pheatmap(all_rpkm,treeheight_row=0,treeheight_col=0,cellwidth = 40,
           display_number =T,cluster_rows=T, cluster_cols=T,
           main=title,
           breaks=breaks,
           colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(length(breaks)-1))
    dev.off()
}


# plot expression for 8 gene panels
plot_all_panels = function(rpkms,gtex_rpkm)
{
    breaks = c(0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,
               5000,6000,7000,8000,9000,10000,15000)
    plot_panel(congenital_myopathy, rpkms, gtex_rpkm, "1_congenital_myopaties.png","Congenital myopaties RPKM",breaks)
  
    breaks = c(0,5,10,20,30,40,50,100)
    plot_panel(congenital_myastenic_syndromes, rpkms, gtex_rpkm, "2_congenital_myastenic_syndromes.png","Congenital myastenic syndromes RPKM",breaks)
  
    breaks = c(0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,110)
    plot_panel(channelopathies, rpkms, gtex_rpkm, "3_channelopathies.png","Channelopathies RPKM",breaks)
  
    breaks = c(0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,2000)
    plot_panel(vacuolar_and_others, rpkms, gtex_rpkm, "4_vacuolar_and_others.png","Vacuolar and others RPKM",breaks)
  
    breaks = c(0,5,10,50,100,500,1000,2000,5000,5200)
    plot_panel(distal_myopathies, rpkms, gtex_rpkm, "5_distal_myopathies.png","Distal myopathies RPKM",breaks)
  
    breaks = c(0,5,10,50,100,500,1000,2000,5000,10000,15000)
    plot_panel(congenital_muscular_dystrophies, rpkms, gtex_rpkm, "6_congenital_muscular_dystrophies.png","Congenital MD RPKM",breaks)
  
    breaks = c(0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,
               5000)
    plot_panel(limb_girdle, rpkms, gtex_rpkm, "7_Limb_girdle_dystrophies.png","Limb_girdle_dystrophies RPKM",breaks)
  
    breaks = c(0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,2000)
    plot_panel(muscular_dystrophies, rpkms, gtex_rpkm, "8_Muscular_dystrophies.png","Muscular_dystrophies RPKM",breaks)
}

expression_fibroblasts = function()
{
    setwd("~/Desktop/project_muscular/Fibroblasts/")
    fibroblast5 = load_rpkm_counts("fibroblast5.rpkm")
    fibroblast6 = load_rpkm_counts("fibroblast6.rpkm")
    fibroblast7 = load_rpkm_counts("fibroblast7.rpkm")
    fibroblast8 = load_rpkm_counts("fibroblast8.rpkm")
    
    fibroblast6$external_gene_name=NULL
    fibroblast7$external_gene_name=NULL
    fibroblast8$external_gene_name=NULL
    
    fibroblasts = merge_row_names(fibroblast5,fibroblast6)
    fibroblasts = merge_row_names(fibroblasts,fibroblast7)
    fibroblasts = merge_row_names(fibroblasts,fibroblast8)
    
    breaks = c(0,5,10,50,100,500,1000)
    
    plot_panel(congenital_muscular_dystrophies,fibroblast8,gtex_rpkm,"fibroblast8.congenital_m_dystrophies.png",
               "Fibroblast8 expr(rpkm), congenital md panel",breaks)
    
    plot_panel(congenital_muscular_dystrophies,fibroblasts,gtex_rpkm,"fibroblasts.congenital_m_dystrophies.png",
               "Fibroblasts expr(rpkm), congenital md panel",breaks)
}

expression_rpkm_muscle2 = function()
{
    setwd("~/Desktop/project_muscular/Muscle2/expression/")
    s62_AF_S5 = load_rpkm_counts("62_AF_S5.rpkm")
    write.table(s62_AF_S5,"62_AF_S5.rpkms.txt",quote=F,sep = "\t")
    
    s1258_AC_A79 = load_rpkm_counts("1258-AC-A79.rpkm")
    s1275_BK_B225 = load_rpkm_counts("1275-BK-B225.rpkm")
    s1388_MJ_M219 = load_rpkm_counts("1388-MJ-M219.rpkm")
    
    s1258_AC_A79$external_gene_name = NULL
    s1275_BK_B225$external_gene_name = NULL
    s1388_MJ_M219$external_gene_name = NULL
    
    s62_AF_S5 = merge_row_names(s62_AF_S5,s1258_AC_A79)
    s62_AF_S5 = merge_row_names(s62_AF_S5,s1275_BK_B225)
    s62_AF_S5 = merge_row_names(s62_AF_S5,s1388_MJ_M219)
    
    write.table(s62_AF_S5,"62_AF_S5.rpkms.controls.txt",quote=F,sep = "\t")
    #sometimes contains duplicate entries, delete them
    s62_AF_S5 = read.delim("62_AF_S5.rpkms.controls.txt", row.names=1, stringsAsFactors=F)
    
    plot_all_panels(s62_AF_S5,gtex_rpkm)
}

expression_rpkm_sample5 = function()
{
    
    setwd("~/Desktop/project_muscular/counts/muscular_filtered/")
    
    #all samples but 1
    samples = read.table("samples.txt", quote="\"", comment.char="", stringsAsFactors=F)
    for (sample in unlist(samples))
    {
        temp = read.delim(paste0(sample,".rpkm.counts.txt"), quote="\"",stringsAsFactors=F, row.names=1)
        temp$Chr=NULL
        temp$Start=NULL
        temp$End=NULL
        temp$Strand=NULL
        temp$Length=NULL
        
        counts = merge(counts,temp,by.x="row.names",by.y="row.names")
        row.names(counts)=counts$Row.names
        counts$Row.names=NULL
    }
    samples = c("Muscle5","Fibroblast5","Myotubes5")
    #counts = counts[c(samples,"symbol")]
    
    ensembl_w_description <- read.delim2("~/Desktop/reference_tables/ensembl_w_description.txt", row.names=1, stringsAsFactors=FALSE)
    counts = merge(counts,ensembl_w_description,by.x="row.names",by.y="row.names")
    row.names(counts)=counts$Row.names
    counts$Row.names=NULL
    
    counts=subset(counts,Muscle5 !=0 & Fibroblast5 !=0 & Myotubes5!=0)
    counts$Gene_description = NULL
    
    write.table(counts,"sample5.txt",quote=F,sep = "\t")
    
    sample5 <- read.csv("sample5.txt", row.names=1, sep="", stringsAsFactors=F)
    
    x=sample5[samples]
    group = factor(c(1,1,1))
    y=DGEList(counts=x,group=group,genes=row.names(x),remove.zeros = F)
    plotMDS(y)
    #generate counts in the other way
    rpkms = rpkm(y,sample5$Length)
    
    rpkms = merge(rpkms,ensembl_w_description,by.x="row.names",by.y="row.names")
    row.names(rpkms)=rpkms$Row.names
    rpkms$Row.names=NULL
    rpkms$Gene_description=NULL
    
    write.table(rpkms,"sample5.rpkms.txt",quote=F,sep = "\t")
    
    rpkms <- read.delim("~/Desktop/project_muscular/counts/muscular_filtered/sample5.rpkms.txt", stringsAsFactors=FALSE)
    
    MuscleGeneRPKM <- read.csv("~/Desktop/project_muscular/counts/muscular_filtered/MuscleGeneRPKM.txt", sep="", stringsAsFactors = F)
    
    breaks = c(0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,500,1000,2000,3000,4000)
    plot_panel(congenital_myopathy, rpkms, MuscleGeneRPKM, "congenital_myopaties",breaks)
    
    breaks = c(0,5,10,50,100,200)
    plot_panel(congenital_myastenic_syndromes, rpkms, MuscleGeneRPKM, "congenital_myastenic_syndromes",breaks)
    
    breaks = c(0,5,10,50)
    plot_panel(channelopathies, rpkms, MuscleGeneRPKM, "channelopathies",breaks)
    
    breaks = c(0,5,10,50,100,500,1000,1500)
    plot_panel(vacuolar_and_others, rpkms, MuscleGeneRPKM, "vacoular_and_others",breaks)
    
    breaks = c(0,5,10,50,100,500,1000,2000,5000,5200)
    plot_panel(distal_myopathies, rpkms, MuscleGeneRPKM, "distal_myopathies",breaks)
    
    breaks = c(0,5,10,50,100,500,1000,2000,5000,10000,15000)
    plot_panel(congenital_muscular_dystrophies, rpkms, MuscleGeneRPKM, "congenital_muscular_dystrophies",breaks)
    
    breaks = c(0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,500,1000,1400)
    #2000,5000,5200)
    plot_panel(limb_girdle, rpkms, MuscleGeneRPKM, "Limb_girdle_dystrophies",breaks)
  
    breaks = c(0,5,10,50,100,500,1000)
    plot_panel(muscular_dystrophies, rpkms, MuscleGeneRPKM, "Muscular_dystrophies",breaks)
}

fibroblast8 = function()
{
  setwd("~/Desktop/project_muscular/counts/mh_unfiltered/")
  samples = 
  mh = read.delim("mh.txt", stringsAsFactors=F, row.names=1)
  fibroblast8 = read.delim("fibroblast8.counts", stringsAsFactors=F, row.names=1)
  
  counts = merge(mh,fibroblast8, by.x='row.names',by.y='row.names')
  row.names(counts) = counts$Row.names
  counts$Row.names = NULL
  
  congenital_muscular_dystrophies=c("LAMA2", "COL6A1", "COL6A2", "COL6A3", 
                                    "SPEN1", "FHL1", "ITGA7", "DNM2","TCAP", "LMNA", "FKTN", 
                                    "POMT1", "POMT2", "FKRP", "POMGNT1", "ISPD", "GTDC2", "B3GNT1", 
                                    "POMGNT1", "GMPPB", "LARGE", "DPM1", "DPM2", "ALG13", "B3GALNT2", 
                                    "TMEM5",  "POMK", "CHKB", "ACTA1", "TRAPPC11")
  
  
  panel.counts = counts[counts$symbol %in% gene_panel,]
  
  
  row.names(all.rpkm) = all.rpkm$external_gene_name
  all.rpkm$external_gene_name=NULL
  
  library(pheatmap)
  library(RColorBrewer)
  png(paste0(title,".png"),res=100,width=500)
  pheatmap(all.rpkm,treeheight_row=0,treeheight_col=0,cellwidth = 40,
           display_number =T,cluster_rows=T, cluster_cols=T,
           main=paste0(title," RPKM"),
           breaks=breaks,
           colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(length(breaks)-1))
  dev.off()
  
}

init = function()
{
    library(edgeR)
    library(pheatmap)
    library(RColorBrewer)
    source("~/Desktop/bioscripts/rnaseq.load_rpkm_counts.R")
    setwd("~/Desktop/project_muscular/")
    
    gene_lengths = read.delim("~/Desktop/project_muscular/reference/gene_lengths.txt", stringsAsFactors=F, row.names=1)
    gtex_rpkm = read.csv("~/Desktop/project_muscular/reference/gtex.muscle_gene.rpkm", sep="", stringsAsFactors = F)
    
    ensembl_w_description = read.delim2("~/Desktop/reference_tables/ensembl_w_description.txt", row.names=1, stringsAsFactors=F)
    
    
    congenital_myopathy = c("ACTA1","TPM3", "TPM2", "TNNT1", "NEB", "LMOD3", "KBTBD13", "CFL2", 
                            "KLHL40", "KLHL41", "MYO18B", "RYR1", "CACNAS1", "STAC3", "ORAI1", 
                            "SITM1", "SEPN1", "CCDC78", "BIN1", "DNM2", "MTM1", "MTMR14", 
                            "SPEG", "PTPLA", "TTN", "MYH7", "MYH2", "CNTN1", "MEGF10","ZAK")
    
    congenital_muscular_dystrophies=c("LAMA2", "COL6A1", "COL6A2", "COL6A3", 
                                      "SPEN1", "FHL1", "ITGA7", "DNM2","TCAP", "LMNA", "FKTN", 
                                      "POMT1", "POMT2", "FKRP", "POMGNT1", "ISPD", "GTDC2", "B3GNT1", 
                                      "POMGNT1", "GMPPB", "LARGE", "DPM1", "DPM2", "ALG13", "B3GALNT2", 
                                      "TMEM5",  "POMK", "CHKB", "ACTA1", "TRAPPC11")
    
    congenital_myastenic_syndromes = c("CHRNA1", "CHRNB1", "CHRND", "CHRNE", "CHRNA1", "RAPSN", "CHAT", "COLQ", "MUSK", "DOK7", 
                                       "AGRN", "GFPT1", "DAPGT1", "LAMB2", "SCN4A", "CHRNG", "PLEC", "ALG2", "ALG14", "SYT2", "PREPL")
    
    channelopathies = c("DMPK", "ZNF9", "CAV3", "HSPG2", "SERCA1", "CLCN1", "SCN4A", "CACNA1S", "CACNA1A", "KCNE3", "KCNA1", "KCNJ18")
    
    vacuolar_and_others = c("LAMP2", "VMA21", "CLN3", "PABPN1", "TTN", "PEC1", "GDF8", "ACVR1", "CAV3", "FHL1", "ICSU","MCOLN1")
    
    limb_girdle = c("MYOT", "LMNA", "CAV3", "DNAJB6", "DES", "TNPO3", "HNRPDL","CAPN3", "DYSF", 
                    "SGCG", "SGCA", "SGCB", "SGCD", "TCAP", "TRIM32", "FKRP", "TTN", "POMT1", "ANO5",
                    "FKTN", "POMT2", "POMGNT1", "PLEC", "DAG1", "TRAPPC11", "GMPPB", "DPM3", "ISPD", "VCP", "LIMS2", "PLEC1", "GAA")
   
    distal_myopathies=c("DYSF", "TTN", "GNE", "MYH7", "MATR3", "TIA1", "MYOT", "NEB", "LDB3", "ANO5", "KLHL9", "DNM2", "FLNC", 
                        "VCP", "CRYAB", "DES", "SPEN1", "LDB3", "FLNC", "BAG3", "TRIM63", "TRIM54")
    
    muscular_dystrophies= c("DMD", "EMD", "FHL1", "LMNA", "SYNE1", "SYNE2", "TMEM43", "TOR1AIP1","DUX4", "SMCHD1", "PTRF")
}
