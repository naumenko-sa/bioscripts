# https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4387895/
# https://cran.r-project.org/web/packages/gplots/gplots.pdf
# ftp://cran.r-project.org/pub/R/web/packages/pheatmap/pheatmap.pdf
# http://www.stat.purdue.edu/~jrounds/weake/sig_10_40/03_edgeR_examine_time.R
# https://cgrlucb.wikispaces.com/edgeR+fall2013

#plot go pictures
go_analysis = function (lrt,prefix)
{
  go = goana(lrt,species="Hs")
  
  for(on in c("BP","CC","MF"))
  {
    go.up = topGO(go,on=on,sort="Up",n=4)
    go.up$log2pvalue = -log2(go.up$P.Up)
    
    go.down = topGO(go,ont=on,sort="Down",n=4)
    go.down$log2pvalue = -log2(go.down$P.Down)
    
    png(paste0(prefix,".HI.",on,".png"),width=1000)
    barplot(go.up$log2pvalue,horiz=T,
            xlab = "-Log2 (Pvalue)",col = "cornflowerblue",cex.axis=1.5,cex.lab = 1.5)
    labels = go.up$Term
    text(x=rep(0.2,4),y=c(0.65,1.85,3.05,4.25),pos=rep(4,1),labels=labels,cex=1.5,font=2)
    dev.off()
    
    png(paste0(prefix,".LO.",on,".png"),width=1000)
    barplot(go.down$log2pvalue,horiz=T,
            xlab = "-Log2 (Pvalue)",col = "cornflowerblue",cex.axis=1.5,cex.lab = 1.5)
    labels = go.down$Term
    text(x=rep(0.2,4),y=c(0.65,1.85,3.05,4.25),pos=rep(4,1),labels=labels,cex=1.5,font=2)
    dev.off()
  }
}

#pathway analysis
kegg_analysis = function (lrt,prefix)
{
  kegg = kegga(lrt,species="Hs")
  kegg.up = topKEGG(kegg,sort="Up",number = 4)
  kegg.up$log2pvalue = -log2(kegg.up$P.Up)
  kegg.down = topKEGG(kegg,sort="Down",number = 4)
  kegg.down$log2pvalue = -log2(kegg.down$P.Down)

  png(paste0(prefix,".HI.kegg.png"),width=1000)
  barplot(kegg.up$log2pvalue,horiz=T,
        xlab = "-Log2 (Pvalue)",col = "cornflowerblue",cex.axis=1.5,cex.lab = 1.5)
  labels = kegg.up$Pathway
  text(x=rep(0.2,4),y=c(0.65,1.85,3.05,4.25),pos=rep(4,1),labels=labels,cex=1.5,font=2)
  dev.off()

  png(paste0(prefix,".LO.kegg.png"),width=1000)
  barplot(kegg.down$log2pvalue,horiz=T,
        xlab = "-Log2 (Pvalue)",col = "cornflowerblue",cex.axis=1.5,cex.lab = 1.5)
  labels = kegg.down$Pathway
  text(x=rep(0.2,4),y=c(0.65,1.85,3.05,4.25),pos=rep(4,1),labels=labels,cex=1.5,font=2)
  dev.off()
}

calc_de = function(all_counts,samples,prefix,filter)
{
    ####    test: ###################
    #samples = c("SG511_ven_hi_2_26","SG511_ven_hi_4_13","SG511_ven_hi_4_27",
    #            "SG523_ven_hi_2_27","SG523_ven_hi_4_10","SG523_ven_hi_4_24",
    #            "SG511_ven_lo_2_26","SG511_ven_lo_4_13","SG511_ven_lo_4_27",
    #            "SG523_ven_lo_2_27","SG523_ven_lo_4_10","SG523_ven_lo_4_24");
    #################################
    
    #samples = c("SG523_ven_lo_2_27","SG523_ven_lo_4_10","SG523_ven_lo_4_24",
    #           "SG523_ven_hi_2_27","SG523_ven_hi_4_10","SG523_ven_hi_4_24")
  
    #all_counts = counts
    #samples = c("G432","G511", "G472","G523", 
    #            "G440","G481", "G510","G564")
    
    
    #samples=samples.523.3points
    #prefix="523.3points_nob_cpm1"
    #filter=1
    
    n_samples = length(samples)
    group=factor(c(rep(1,n_samples/2),rep(2,n_samples/2)))
    #patient = factor(c("511","511","511","523","523","523",
    #                  "511","511","511","523","523","523"))
    
    x=all_counts[samples]
    y=DGEList(counts=x,group=group,genes=row.names(x),remove.zeros = T)

    plotMDS(y)
    #filter - 1 or 0.5
    #filter=1
    keep=rowSums(cpm(y)>filter) >= n_samples/2
    y=y[keep,,keep.lib.sizes=F]

    #necessary for goana
    idfound = y$genes$genes %in% mappedRkeys(org.Hs.egENSEMBL)
    y = y[idfound,]
    
    egENSEMBL=toTable(org.Hs.egENSEMBL)
    m = match (y$genes$genes,egENSEMBL$ensembl_id)
    y$genes$EntrezGene = egENSEMBL$gene_id[m]
    egSYMBOL = toTable(org.Hs.egSYMBOL)
    m = match (y$genes$EntrezGene,egSYMBOL$gene_id)
    y$genes$Symbol = egSYMBOL$symbol[m]
    
    #remove duplications - just 1 gene in this dataset
    o = order(rowSums(y$counts),decreasing = T)
    y = y[o,]
    d = duplicated(y$genes$Symbol)
    dy = y[d,]$genes
    y = y[!d,]
    nrow(y)
    
    y$samples$lib.size = colSums(y$counts)
    rownames(y$counts) = y$genes$EntrezGene 
    rownames(y$genes) = y$genes$EntrezGene
    y$genes$EntrezGene = NULL
    
    #normalization for RNA composition (2.7.3)
    y=calcNormFactors(y)

#nc=cpm(y,normalized.lib.sizes=F)
#write.table(nc,"filtered.normalized_counts.txt",col.names=NA)

    png(paste0(prefix,".png"))
    plotMDS(y,las=1)
    dev.off()

    design=model.matrix(~group)

    y=estimateDisp(y,design)
    
    fit=glmFit(y,design)
    lrt=glmLRT(fit)

    #o=order(lrt$table$PValue)
    #cpm(y)[o[1:10],]
    #write.table(cpm(y)[o[1:12566],],"allgenes.cpm")
    
    efilename="all_hi_vs_lo.genes.txt"
    write.table(topTags(lrt,p.value=0.05,n=10000),efilename,quote=F)

    de_results = read.csv(efilename, sep="", stringsAsFactors=FALSE)
    s_rownames = row.names(de_results)
    #setnames(de_results,"genes","ensembl_gene_id")
    #de_results = lrt$table
    
    gene_descriptions = read.delim2(paste0(reference_tables_path,"/ensembl_w_description.txt"), stringsAsFactors=FALSE)
    
    de_results = merge(de_results,gene_descriptions,by.x="genes",by.y="ensembl_gene_id",all.x=T)
    #de_results = rename(de_results,c("Row.names"="ensembl_gene_id"))
    de_results = merge(de_results,x,by.x = "genes", by.y="row.names",all.x=T)
    de_results = de_results[order(de_results$PValue),]
    rownames(de_results) = s_rownames
    
    result_file=paste0(prefix,".txt")
    write.table(de_results,result_file,quote=F,sep=';')

    plot_heatmap_separate(all_counts,samples,de_results,prefix)    
    #return(de_results)
    
    go_analysis(lrt,prefix)
    
    #kegg_analysis(lrt,prefix)
    
}
 
plot_heatmap_separate = function(counts,samples,de_results,prefix)
{
  logcpm = cpm(counts,prior.count=1,log=T)
  #cpm0  = log(cpm(y$counts+1))
  top_genes_cpm = logcpm[de_results$genes,]
  top_genes_cpm = top_genes_cpm[,samples]
  rownames(top_genes_cpm) = de_results$external_gene_name
  
  png(paste0(prefix,".heatmap.png"))
  
  pheatmap(top_genes_cpm,scale="row",treeheight_row=0,treeheight_col=0,cellwidth = 20)
  
  dev.off()
}   

plot_heatmap = function (counts,de_results,prefix)
{
    logcpm = cpm(counts,prior.count=1,log=T)
    #cpm0  = log(cpm(y$counts+1))
    top_genes_cpm = logcpm[de_results$genes,]
    rownames(top_genes_cpm) = de_results$external_gene_name
    
    png(paste0(prefix,".heatmap.png"),width=2000)
    
    pheatmap(top_genes_cpm,scale="row",treeheight_row=0)
    
    ph=pheatmap(t(top_genes_cpm), kmeans_k=2, scale="column", show_rownames=F,
             treeheight_row = 0, treeheight_col = 0, fontsize = 20,
             cellheight = 60, cellwidth = 20,rot=90,cluster_rows = T,
             show_colnames = F)
    
    #ph=pheatmap(t(top_genes_cpm)[c(3,4,1,2),],kmeans_k = 2,show_rownames = F,treeheight_row = 0,
    #         treeheight_col = 0, cellheight = 60, cellwidth = 20,
    #         fontsize=20, scale="column", cluster_rows = T,
    #         color = colorRampPalette(c("green", "black", "red"))(20),
    #         show_colnames = F, rot=90)
    labels = ph$tree_col$labels[ph$tree_col$order]
    #start = 0.88 - length(labels)*0.01
             #filename="SG523.heatmap.png")
    start = 0.555 - length(labels)*0.01
    grid.text(labels,
              seq(start,start+(length(labels)-1)*0.01,0.01),
              rep(0.638,length(labels)),
              rot=rep(90,length(labels)),
              just = "left")
    
    dev.off()
}   

#using the section 4.2 of edgeR manual
calc_de_w_batch_effect = function(counts,samples,prefix,filter)
{
    #samples = c("SG511_ven_lo_4_13","SG511_ven_lo_4_27",
    #          "SG523_ven_lo_2_27","SG523_ven_lo_4_10","SG523_ven_lo_4_24",
    #          "SG511_ven_hi_4_13","SG511_ven_hi_4_27",
    #          "SG523_ven_hi_2_27","SG523_ven_hi_4_10","SG523_ven_hi_4_24"
    #          );
    #counts=read.delim("combined.counts",row.names="id")
    #prefix = "test"
    #samples = samples.511.2
    #counts = all_counts
    #prefix = "523.3points_wb_cpm1"
    #samples = samples.523.3points
    #filter=1
  
    x = counts[samples]
    n_samples = length(samples)
    
    treat = factor(substring(colnames(x),11,12),levels=c("lo","hi"))
    
    #try 3 times
    time=factor(substring(colnames(x),14,17))
    
    patient = factor(substring(colnames(x),3,5))
    
    #y=DGEList(counts=x,group=treat)
    y=DGEList(counts=x,group=treat,genes=row.names(x),remove.zeros = T)
    
    #>=0.5 filter =1 or 0.5calc_de_w_batch_effect(all_counts,samples.523.3points,"523.3points_wb_cpm0.5",0.5)calc_de_w_batch_effect(all_counts,samples.523.3points,"523.3points_wb_cpm0.5",0.5)calc_de_w_batch_effect(all_counts,samples.523.3points,"523.3points_wb_cpm0.5",0.5)
    keep = rowSums(cpm(y)>filter) >= n_samples/2
    table(keep)
 
    y = y [keep, ,keep.lib.sizes=F]   
    
    #necessary for goana
    idfound = y$genes$genes %in% mappedRkeys(org.Hs.egENSEMBL)
    y = y[idfound,]
    
    egENSEMBL=toTable(org.Hs.egENSEMBL)
    m = match (y$genes$genes,egENSEMBL$ensembl_id)
    y$genes$EntrezGene = egENSEMBL$gene_id[m]
    egSYMBOL = toTable(org.Hs.egSYMBOL)
    m = match (y$genes$EntrezGene,egSYMBOL$gene_id)
    y$genes$Symbol = egSYMBOL$symbol[m]
    
    #remove duplications - just 1 gene in this datasetcalc_de_w_batch_effect(all_counts,samples.523.3points,"523.3points_wb_cpm0.5",0.5)
    o = order(rowSums(y$counts),decreasing = T)
    y = y[o,]
    d = duplicated(y$genes$Symbol)
    dy = y[d,]$genes
    y = y[!d,]
    nrow(y)
    
    y$samples$lib.size = colSums(y$counts)
    rownames(y$counts) = y$genes$EntrezGene 
    rownames(y$genes) = y$genes$EntrezGene
    y$genes$EntrezGene = NULL
    
    y = calcNormFactors(y)
    
    y$samples
    
    plotMDS(y)
    
    #design = model.matrix(~time+time:treat)
    #logFC = predFC(y,design,prior.count=1,dispersion=0.05)
    
    #cor(logFC[,6:10])    
    
    design = model.matrix(~time+treat)
    rownames(design) = colnames(y)
    design
    
    y = estimateDisp(y,design,robust=T)
    
    y$common.dispersion
    
    plotBCV(y)
    
    fit=glmQLFit(y,design,robust=T)
    plotQLDisp(fit)
    
    #check DE for time - 623 genes
    #qlf=glmQLFTest(fit,coef=2:3)
    #topTags(qlf)
    #FDR = p.adjust(qlf$table$PValue, method="BH")
    #sum(FDR<0.05)
    #go_analysis(qlf,"5points.time")
    
    #de for HI/Lo
    qlf = glmQLFTest(fit)
    efilename="tmp.txt"
    write.table(topTags(qlf,p.value=0.05,n=1000),efilename,quote=F)
    de_results = read.csv(efilename, sep="", stringsAsFactors=FALSE)

    gene_descriptions = read.delim2(paste0(reference_tables_path,"/ensembl_w_description.txt"), stringsAsFactors=FALSE)
    
    de_results = merge(de_results,gene_descriptions,by.x="genes",by.y="ensembl_gene_id",all.x=T)
    #de_results = rename(de_results,c("Row.names"="ensembl_gene_id"))
    de_results = merge(de_results,x,by.x = "genes", by.y="row.names",all.x=T)
    
    write.table(de_results,paste0(prefix,".txt"),quote=F,sep=';')
    
    plot_heatmap_separate(all_counts,samples,de_results,prefix)
    
    go_analysis(qlf,prefix)
    #kegg_analysis(qlf,paste0(prefix,".batch"))
    
    #dt = decideTestsDGE(qlf)
    #summary(dt)
    
    #top=rownames(topTags(qlf))
    #View(cpm(y)[top,])
}

init = function()
{
  library(edgeR)
  library(stringr)
  #library(plyr)
  library(GO.db)
  library(org.Hs.eg.db)
  library(data.table)
  library(pheatmap)
  library(grid)  
  library("VennDiagram")
  reference_tables_path = "~/Desktop/reference_tables"
  setwd("~/Desktop/project_katie_csc")

#samples
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

}

test8samples = function()
{
    counts = read.delim("wG432_DMSO.txt", stringsAsFactors=F, row.names=1)
    #all samples but 1
    samples = read.table("samples.txt", quote="\"", comment.char="", stringsAsFactors=F)
    for (sample in unlist(samples))
    {
        temp = read.delim(paste0(sample,".txt"), quote="\"", 
                          comment.char="",stringsAsFactors=F, row.names=1)
        counts = merge(counts,temp,by.x="row.names",by.y="row.names")
        row.names(counts)=counts$Row.names
        counts$Row.names=NULL
    }
    
    samples = c("G432","G511", "G472","G523", 
                "G440","G481", "G510","G564")
    
    samples = c("G511", "G472","G523", 
                "G440", "G510","G564")
    prefix = "DMSO_6lines"
    all_counts=read.csv("DMSO_counts.txt", row.names=1, sep="", stringsAsFactors=F)
    calc_de(all_counts,samples,prefix,1)
}

samples5 = c("SG511_ven_lo_4_13","SG511_ven_lo_4_27",
          "SG523_ven_lo_2_27","SG523_ven_lo_4_10","SG523_ven_lo_4_24",
          "SG511_ven_hi_4_13","SG511_ven_hi_4_27",
          "SG523_ven_hi_2_27","SG523_ven_hi_4_10","SG523_ven_hi_4_24"
          )

#511.2points
samples.511.2points = c("SG511_ven_lo_4_13","SG511_ven_lo_4_27",
            "SG511_ven_hi_4_13","SG511_ven_hi_4_27")


#523.2points
samples.523.2points = c("SG523_ven_lo_2_27","SG523_ven_lo_4_10",
            "SG523_ven_hi_2_27","SG523_ven_hi_4_10")

#final set - 4 points from two samples:
samples.523.2_511.2 = c("SG523_ven_lo_2_27","SG523_ven_lo_4_10",
                        "SG511_ven_lo_4_13","SG511_ven_lo_4_27",
                        "SG523_ven_hi_2_27","SG523_ven_hi_4_10",
                        "SG511_ven_hi_4_13","SG511_ven_hi_4_27")

samples.523.3points = c("SG523_ven_lo_2_27","SG523_ven_lo_4_10","SG523_ven_lo_4_24",
                        "SG523_ven_hi_2_27","SG523_ven_hi_4_10","SG523_ven_hi_4_24")

calc_de_w_batch_effect(all_counts,samples5,"5points_wb_cpm1",1)

calc_de_w_batch_effect(all_counts,samples.523.2,"523.2points")
calc_de_w_batch_effect(all_counts,samples.523.2_511.2,"4points_wb_cpm1",1)
calc_de_w_batch_effect(all_counts,samples.523.3points,"523.3points")

calc_de_w_batch_effect(all_counts,samples.523.3points,"523.3points_wb_cpm0.5",0.5)
calc_de_w_batch_effect(all_counts,samples.523.3points,"523.3points_wb_cpm1",1)

calc_de_w_batch_effect(all_counts,samples.523.3points,"523.3points_wb_cpm2",2)
calc_de_w_batch_effect(all_counts,samples.523.2points,"523.2points_wb_1cpm",1)
calc_de_w_batch_effect(all_counts,samples.523.2points,"523.2points_wb_0.5cpm",0.5)

calc_de_w_batch_effect(all_counts,samples.511.2points,"511.2points_wb_cpm0.5",0.5)


calc_de(all_counts,samples.523.3points,"523.3points_nob_cpm1",1)
calc_de(all_counts,samples.523.3points,"523.3points_nob_cpm0.5",0.5)
calc_de(all_counts,samples.523.2points,"523.2points_nob_cpm0.5",0.5)

r523.2_511.2points= read.csv("SG523/3points/523.2_511.2_w_batch_effect_correction.txt",header=T,sep=";")
r523.2points = read.csv("~/Dropbox/project_katie_csc/523.2points_w_batch_effect_correction.txt", header=T, sep=";")

r523.3points = read.csv("de_results/SG523/3points/wb_cpm0.5/523.3points_wb_cpm0.5.txt", header=T, sep=";")
r523.2points_no_batch = read.csv("~/Dropbox/project_katie_csc/523.2points.txt", header=T, sep=";")

r523.6points = read.csv("~/Dropbox/project_katie_csc/all_sample_w_batch_effect1.txt", header=T, sep=";")
r5points = read.csv("~/Dropbox/project_katie_csc/5points_w_batch_effect.txt", header=T, sep=";")
r511.2points = read.csv("de_results/SG511/2points/wb_0.5cpm/511.2points_wb_cpm0.5.txt",header=T,sep=";")

overlap = calculate.overlap(x=list("511.2points" = r511.2points$Symbol,
                                    "523.3points" = r523.3points$Symbol))

overlap = calculate.overlap(x=list("511.2points_batch" = r511.2points$Symbol,
                                   "5points_batch" = r5points$Symbol,
                                   "523.2points" = r523.2points$Symbol))

overlap = calculate.overlap(x=list("523.3points" = r523.3points$Symbol,
                                   "523.2points" = r523.2points_no_batch$Symbol))

overlap = calculate.overlap(x=list("523.2points" = r523.2points$Symbol,
                                   "511.2points" = r511.2points$Symbol))


png("511.2points_523.3points.png",width=1200)
venn.plot=draw.pairwise.venn(length(overlap$a1),
                             length(overlap$a2),
                             length(overlap$a3),
                             c("SG511_2points","SG523_3points"),fill=c("blue","red"),
                             lty="blank",
                             cex = 4, cat.cex=4,cat.pos = 0,lwd=c(2,2),
                             ext.length = 0.3,  
                             ext.line.lwd=2,
                             ext.text = F,scaled = F)
dev.off()

png("5points.511_2points.523_2points.png",width=1000)
venn.plot=draw.triple.venn(area1 = length(r5points$Symbol),
                           area2 = length(r523.2points$Symbol),
                           area3 = length(r511.2points$Symbol),
                           n12 = length(intersect(r5points$Symbol,r523.2points$Symbol)),
                           n23 = length(intersect(r523.2points$Symbol,r511.2points$Symbol)),
                           n13 = length(intersect(r5points$Symbol,r511.2points$Symbol)),
                           n123 = length(intersect(intersect(r511.2points$Symbol,r5points$Symbol),r523.2points$Symbol)),
                           category = c("5points","523.2points","511.2points"),fill=c("blue","red","yellow"),
                             lty="blank",
                             cex = 2, cat.cex=2, 
                             ext.length = 0.3,  
                             ext.line.lwd=2,
                             ext.text = F,
                           scaled=F,
                           euler.d=T)
dev.off()

