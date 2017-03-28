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
  
    #test:
    all_counts = counts
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
    max_genes = nrow(x)
    
    logcpm = cpm(all_counts,prior.count=1,log=T)
    t_cpm = cpm(all_counts,prior.count=1,log=F)
    logcpm = logcpm[,samples]
    t_cpm = t_cpm[,samples]
    
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
    #first order by counts to remove duplicated names with 0 counts
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

    png(paste0(prefix,".mds.png"),res = 300,width=2000,height=2000)
    plotMDS(y,las=1)
    dev.off()

    design=model.matrix(~group)

    y=estimateDisp(y,design)
    
    fit=glmFit(y,design)
    lrt=glmLRT(fit)

    efilename=paste0(prefix,".de_genes.txt")
    de_results = topTags(lrt,p.value=0.05,n=max_genes,sort.by="logFC")
    write.table(de_results,efilename,quote=F,row.names=F)

    de_results = read.csv(efilename, sep="", stringsAsFactors=FALSE)
    s_rownames = row.names(de_results)
    #setnames(de_results,"genes","ensembl_gene_id")
    #de_results = lrt$table
    
    gene_descriptions = read.delim2(paste0(reference_tables_path,"/ensembl_w_description.txt"), stringsAsFactors=FALSE)
    
    de_results = merge(de_results,gene_descriptions,by.x="genes",by.y="ensembl_gene_id",all.x=T)
    #de_results = rename(de_results,c("Row.names"="ensembl_gene_id"))
    de_results = merge(de_results,x,by.x = "genes", by.y="row.names",all.x=T)
    #rownames(de_results) = s_rownames
    
    top_genes_cpm = logcpm[de_results$genes,]
    colnames(top_genes_cpm)=paste0(colnames(top_genes_cpm),".log2cpm")
    
    de_results = merge(de_results,top_genes_cpm,by.x = "genes", by.y="row.names",all.x=T)
    de_results = de_results[order(abs(de_results$logFC),decreasing = T),]
    
    colnames(de_results)[1]="Ensembl_gene_id"
    colnames(de_results)[2]="Gene_name"
    de_results$external_gene_name = NULL
    result_file=paste0(prefix,".txt")
    write.table(de_results,result_file,quote=T,row.names=F)
    
    prepare_file_4gsea(all_counts,samples,prefix,gene_descriptions)
    
    plot_heatmap_separate (all_counts,samples,de_results,prefix)
    plot_heatmap_separate (all_counts,samples,de_results,paste0(prefix,".top50genes"),50)
    #return(de_results)
    
    #go_analysis(lrt,prefix)
    
    #kegg_analysis(lrt,prefix)
    
}

prepare_file_4gsea = function(counts,samples,prefix,gene_descriptions)
{
  #prepare a file for GSEA - positive and negative in a separate file
  #http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#Expression_Data_Formats
  #for GSEA it is important to report all genes - genome wide
  
    t_cpm = cpm(counts,prior.count=1,log=F)
    t_cpm = t_cpm[,samples]
  
    result_file=paste0(prefix,".4gsea.txt")
    t_cpm =  merge(t_cpm,gene_descriptions,by.x="row.names",by.y="ensembl_gene_id",all.x=T)
    colnames(t_cpm)[1]="Ensembl_gene_id"
    t_cpm = t_cpm[c("external_gene_name","Ensembl_gene_id",paste0(samples))]
    colnames(t_cpm)[1:2]=c("NAME","DESCRIPTION")
  
    o = order(rowSums(t_cpm[,c(samples)]),decreasing = T)
    t_cpm = t_cpm[o,]
    d = duplicated(t_cpm$NAME)
    dy = t_cpm[d,]$NAME
    t_cpm = t_cpm[!d,]
    nrow(t_cpm)
  
    write.table(t_cpm,result_file,quote=F,row.names = F,sep = "\t")
}

# usually heatmap is a part of a panel - we don't need a title
# rows are not clustered to save alpabetical gene order
# cols are not clustered to save sample order
plot_heatmap = function(prefix,expression_table)
{
    rows = nrow(expression_table)
    cellheight = 10
    res = 300
    filename = paste0(prefix,".",res,"ppi.png")
    
    png(filename,res = res,height=rows * cellheight * 4.5,width=1500)
    pheatmap(expression_table,scale="row",treeheight_row=0,treeheight_col=0,
             display_numbers = T,cellheight = cellheight,cellwidth = 30,
             cluster_rows = F, cluster_cols = F)
    dev.off()
}

# plot separate heatmaps for upregulated and downregulated genes 
# parameters:
# counts - initial row count for all genes
# samples
# de_results
# prefix  
# ntop - how many top genes to plot
plot_heatmap_separate = function(counts,samples,de_results,prefix,ntop = NULL)
{
    logcpm = cpm(counts,prior.count=1,log=T)
    #cpm0  = log(cpm(y$counts+1))
    top_genes_cpm = logcpm[de_results$Ensembl_gene_id,]
    top_genes_cpm = top_genes_cpm[,samples]
    rownames(top_genes_cpm) = de_results$Gene_name
    
    #expressed higher in WNT-dependent cells.
    upregulated_genes = de_results[de_results$logFC<0,]$Gene_name
    downregulated_genes = de_results[de_results$logFC>0,]$Gene_name
    
    if (!is.null(ntop))
    {
        upregulated_genes = head(upregulated_genes,ntop)
        downregulated_genes = head(downregulated_genes,ntop)
    }
    
    #sort genes alphabetically - it is much easier to read heatmap
    upregulated_genes = sort(upregulated_genes)
    downregulated_genes = sort (downregulated_genes)
    
    
    plot_heatmap(paste0(prefix,".left"),top_genes_cpm[upregulated_genes,])
    plot_heatmap(paste0(prefix,".right"),top_genes_cpm[downregulated_genes,])
}   

plot_heatmap1 = function (counts,de_results,prefix)
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

# merge pipeline output counts into 1 count file 
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
    
}

# figure 2C - use lgk samples only, because lgk samples = dmso, same expression
# G481 = G361
figure2C_8samples = function()
{
    setwd("~/Desktop/project_katie_csc/Figure2C/")
    counts = read.csv("LGK_counts.txt", row.names=1, sep="", stringsAsFactors=F)
    prefix = "Figure2C.supp.8cell_lines.log2cpm"
    #columns are clustered like on the left picture G361 - outlier, G432 is close to G564
    samples = c("G472","G511","G523","G432","G564","G440","G510","G361")
    filter=1
    calc_de(counts_lgk,samples,prefix,filter)
}

figure2C_6samples = function()
{
    setwd("~/Desktop/project_katie_csc/Figure2C/")
    counts = read.csv("LGK_counts.txt", row.names=1, sep="", stringsAsFactors=F)
    prefix = "Figure2C.supp.6cell_lines.log2cpm"
    samples = c("G472","G511","G523", 
                "G564","G440","G510")
    filter=1
    calc_de(counts,samples,prefix,1)
}

test_lgk_dmso = function()
{
    setwd("~/Desktop/project_katie_csc/DMSO_LGK_expression/")
    counts_lgk = read.csv("LGK_counts.txt", row.names=1, sep="", stringsAsFactors=F)
    counts_dmso = read.csv("DMSO_counts.txt", row.names=1, sep="", stringsAsFactors=F)
    counts = merge(counts_lgk,counts_dmso,by.x="row.names",by.y="row.names")
    
    row.names(counts)=counts[,1]
    counts$Row.names = NULL
    write.table(counts,"counts.txt",quote=F)
    
    samples = c("G432_lgk","G511_lgk", "G472_lgk","G523_lgk", 
                "G440_lgk","G481_lgk", "G510_lgk","G564_lgk",
                "G432_dmso","G511_dmso", "G472_dmso","G523_dmso", 
                "G440_dmso","G481_dmso", "G510_dmso","G564_dmso")
    
    prefix = "8lines"
    
    calc_de(counts,samples,prefix,1)
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

