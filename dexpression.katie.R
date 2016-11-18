# https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4387895/
# https://cran.r-project.org/web/packages/gplots/gplots.pdf
# ftp://cran.r-project.org/pub/R/web/packages/pheatmap/pheatmap.pdf

calc_de = function(all_counts,samples,result_file)
{
    ####    test: ###################
    #samples = c("SG511_ven_hi_2_26","SG511_ven_hi_4_13","SG511_ven_hi_4_27",
    #            "SG523_ven_hi_2_27","SG523_ven_hi_4_10","SG523_ven_hi_4_24",
    #            "SG511_ven_lo_2_26","SG511_ven_lo_4_13","SG511_ven_lo_4_27",
    #            "SG523_ven_lo_2_27","SG523_ven_lo_4_10","SG523_ven_lo_4_24");
    #################################
    
    #samples = c("SG523_ven_hi_2_27","SG523_ven_hi_4_10","SG523_ven_hi_4_24",
    #          "SG523_ven_lo_2_27","SG523_ven_lo_4_10","SG523_ven_lo_4_24");
    
    #samples = c("SG523_ven_hi_2_27","SG523_ven_hi_4_24",
    #            "SG523_ven_lo_2_27","SG523_ven_lo_4_24");
    
    samples = c("SG523_ven_lo_2_27","SG523_ven_lo_4_10",
                "SG523_ven_hi_2_27","SG523_ven_hi_4_10")

    #samples = c("SG511_ven_lo_4_13","SG511_ven_hi_4_13",
    #            "SG511_ven_lo_2_26","SG511_ven_hi_2_26");

  
    n_samples = length(samples)
    group=factor(c(rep(1,n_samples/2),rep(2,n_samples/2)))
    #patient = factor(c("511","511","511","523","523","523",
    #                  "511","511","511","523","523","523"))
    
    x=all_counts[samples]
    y=DGEList(counts=x,group=group)

    plotMDS(y)
    keep=rowSums(cpm(y)>1) >= n_samples/2
    y=y[keep,,keep.lib.sizes=F]

    #normalization for RNA composition (2.7.3)
    y=calcNormFactors(y)

#nc=cpm(y,normalized.lib.sizes=F)
#write.table(nc,"filtered.normalized_counts.txt",col.names=NA)

    plotMDS(y)

    design=model.matrix(~group)

    y=estimateDisp(y,design)
    
    fit=glmFit(y,design)
    lrt=glmLRT(fit)

    efilename="all_hi_vs_lo.genes.txt"
    write.table(topTags(lrt,p.value=0.05,n=10000),efilename,quote=F)

    de_results = read.csv(efilename, sep="", stringsAsFactors=FALSE)
    #de_results = lrt$table
    
    gene_descriptions = read.delim2(paste0(reference_tables_path,"/ensembl_w_description.txt"), stringsAsFactors=FALSE)
    
    de_results = merge(de_results,gene_descriptions,by.x="row.names",by.y="ensembl_gene_id",all.x=T)
    de_results = rename(de_results,c("Row.names"="ensembl_gene_id"))
    de_results = merge(de_results,x,by.x = "ensembl_gene_id", by.y="row.names",all.x=T)

    result_file="523.2points.txt"
    write.table(de_results,result_file,quote=F,sep=';')
    
    return(de_results)
    
    logcpm = cpm (y,prior.count=2,log=T)
    #cpm0  = log(cpm(y$counts+1))
    s40 = logcpm[de_results$ensembl_gene_id,]
    rownames(s40) = de_results$external_gene_name
    library(pheatmap)
    library(grid)  
    png("SG523_new.png",width=2000)
    ph=pheatmap(t(s40)[c(3,4,1,2),],kmeans_k = 2,show_rownames = F,treeheight_row = 0,
             treeheight_col = 0, cellheight = 60, cellwidth = 20,
             fontsize=20,scale="column", cluster_rows = T,
             color = colorRampPalette(c("green", "black", "red"))(20),
             show_colnames = F, rot=90)
    labels = ph$tree_col$labels[ph$tree_col$order]
             #filename="SG523.heatmap.png")
    grid.text(labels,
              seq(0.09,0.09+(length(labels)-1)*0.01,0.01),
              rep(0.637,length(labels)),
              rot=rep(90,length(labels)),
              just = "left")
    dev.off()
}

calc_de_w_batch_effect = function()
{
    #using the section 4.2 of the manual
    samples = c("SG511_ven_hi_2_26","SG511_ven_hi_4_13","SG511_ven_hi_4_27",
              "SG523_ven_hi_2_27","SG523_ven_hi_4_10","SG523_ven_hi_4_24",
              "SG511_ven_lo_2_26","SG511_ven_lo_4_13","SG511_ven_lo_4_27",
              "SG523_ven_lo_2_27","SG523_ven_lo_4_10","SG523_ven_lo_4_24");

    all_counts=read.delim("combined.counts",row.names="id")
    x = all_counts[samples]
    
    treat = factor(substring(colnames(x),11,12))
    
    time=factor(substring(colnames(x),14,17))
    
    patient = factor(substring(colnames(x),3,5))
    
    y=DGEList(counts=x,group=treat)
    
    keep = rowSums(cpm(y)>1) > 3
    table(keep)
 
    y = y [keep, ,keep.lib.sizes=F]   
    
    y = calcNormFactors(y)
    
    y$samples
    
    plotMDS(y)
    
    design = model.matrix(~time+time:treat)
    logFC = predFC(y,design,prior.count=1,dispersion=0.05)
    
    cor(logFC[,4:6])    
    
    design = model.matrix(~time+treat)
    rownames(design) = colnames(y)
    
    y = estimateDisp(y,design,robust=T)
    
    y$common.dispersion
    
    plotBCV(y)
    
    fit=glmQLFit(y,design,robust=T)
    plotQLDisp(fit)
    
    qlf=glmQLFTest(fit,coef=2:3)
    topTags(qlf)
    
    FDR = p.adjust(qlf$table$PValue, method="BH")
    sum(FDR<0.05)
    
    qlf = glmQLFTest(fit)
    
    result=topTags(qlf)
    
    efilename="all_hi_vs_lo.genes.txt"
    write.table(topTags(qlf,p.value=0.05,n=10000),efilename,quote=F)
    de_results = read.csv(efilename, sep="", stringsAsFactors=FALSE)

    gene_descriptions = read.delim2(paste0(reference_tables_path,"/ensembl_w_description.txt"), stringsAsFactors=FALSE)
    
    de_results = merge(de_results,gene_descriptions,by.x="row.names",by.y="ensembl_gene_id",all.x=T)
    de_results = rename(de_results,c("Row.names"="ensembl_gene_id"))
    de_results = merge(de_results,x,by.x = "ensembl_gene_id", by.y="row.names",all.x=T)
    
    write.table(de_results,"all_sample_w_batch_effect.txt",quote=F,sep=';')
    
    top = row.names(topTags(qlf))    
    View(cpm(y)[top,])
    
    
    
    dt = decideTestsDGE(qlf)
    summary(dt)
    
    isDE= as.logical(dt)
    DEnames = rownames(y)[isDE]
    
    plotSmear(qlf,de.tags=DEnames)
    abline(h=c(-1,1),col="blue")
    
    gene_descriptions = read.delim2(paste0(reference_tables_path,"/ensembl_w_description.txt"), stringsAsFactors=FALSE)
    final_set = gene_descriptions[gene_descriptions$ensembl_gene_id %in% DEnames,]
    final_counts = x[DEnames,]
  }

library(edgeR)
library(stringr)
library(plyr)
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
samples
#SG523_ven_lo_4_10
#SG523_ven_lo_4_24

all_counts=read.delim("combined.counts",row.names="id")

#no outliers
#just 523
samples_all = c("SG511_ven_hi_2_26","SG511_ven_hi_4_13","SG511_ven_hi_4_27",
                "SG523_ven_hi_2_27","SG523_ven_hi_4_10","SG523_ven_hi_4_24",
                "SG511_ven_lo_2_26","SG511_ven_lo_4_13","SG511_ven_lo_4_27",
                "SG523_ven_lo_2_27","SG523_ven_lo_4_10","SG523_ven_lo_4_24");
result_all = calc_de(all_counts,samples_all)

samples523 = c("SG523_ven_hi_2_27","SG523_ven_hi_4_10","SG523_ven_hi_4_24",
            "SG523_ven_lo_2_27","SG523_ven_lo_4_10","SG523_ven_lo_4_24");
results523=calc_de(all_counts,samples523,"523.txt")

#no result
samples511 = c("SG511_ven_hi_2_26","SG511_ven_hi_4_13","SG511_ven_hi_4_27",
            "SG511_ven_lo_2_26","SG511_ven_lo_4_13","SG511_ven_lo_4_27")
samples511 = c("SG511_ven_hi_4_13","SG511_ven_hi_4_27",
               "SG511_ven_lo_4_13","SG511_ven_lo_4_27")
samples511_4_13 = c("SG511_ven_hi_4_13","SG511_ven_lo_4_13","SG511_ven_lo_4_27")
results511 = calc_de(all_counts,samples511,"511.txt")


#http://www.stat.purdue.edu/~jrounds/weake/sig_10_40/03_edgeR_examine_time.R
cpm0 = log(cpm(y$counts+1))
colSums(cpm0)
