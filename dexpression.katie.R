# https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

calc_de = function(all_counts,samples)
{
    ####    test: ###################
    samples = c("SG511_ven_hi_2_26","SG511_ven_hi_4_13","SG511_ven_hi_4_27",
                "SG523_ven_hi_2_27","SG523_ven_hi_4_10","SG523_ven_hi_4_24",
                "SG511_ven_lo_2_26","SG511_ven_lo_4_13","SG511_ven_lo_4_27",
                "SG523_ven_lo_2_27","SG523_ven_lo_4_10","SG523_ven_lo_4_24");
    #################################
    
    samples = c("SG523_ven_hi_2_27","SG523_ven_hi_4_10","SG523_ven_hi_4_24",
                "SG511_ven_hi_4_13",
              "SG523_ven_lo_2_27","SG523_ven_lo_4_10","SG523_ven_lo_4_24",
              "SG511_ven_lo_4_13");

    samples = c("SG511_ven_hi_2_26","SG511_ven_hi_4_13","SG511_ven_hi_4_27",
                "SG511_ven_lo_2_26","SG511_ven_lo_4_13","SG511_ven_lo_4_27");
    
    samples = c("SG511_ven_lo_4_13","SG511_ven_hi_4_13",
                "SG511_ven_lo_2_26","SG511_ven_hi_2_26");
    
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

    return(de_results)
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
    
    keep = rowSums(cpm(y)>1) > 2
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
    
    topTags(qlf)

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
results523=calc_de(all_counts,samples)

samples_test = c("SG511_ven_hi_4_13",
                "SG523_ven_hi_2_27","SG523_ven_hi_4_10",
                "SG511_ven_lo_4_13",
                "SG523_ven_lo_2_27","SG523_ven_lo_4_10")