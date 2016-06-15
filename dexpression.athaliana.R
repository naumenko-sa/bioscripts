# 2016-06-15: DE of miRNA for A.thaliana project

# https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
library(edgeR)
setwd("~/Desktop/project_micro_RNA_athaliana/dexpression")

# 1T4T1 is absent

all_counts=read.delim("counts_mirna.tsv",row.names="mirna")

de_comparison = function(x,group,output)
{
    #test:
    #x=all_counts[c("X12T3R1","X12T3R2","X12T3R3","X12T4R1","X12T4R2","X12T4R3")]
    #group=factor(c(1,1,1,2,2,2))
    #output="test1"
  
    y=DGEList(counts=x,group=group)
    y=calcNormFactors(y)
    design=model.matrix(~group)
    y=estimateDisp(y,design)
    write.table(y$counts,paste0(output,".normalized_counts.txt"),col.names=NA)
    fit=glmFit(y,design)
    lrt=glmLRT(fit)
    #prints top 50 genes with p.value<0.05, check if there are more
    write.table(topTags(lrt,p.value=0.05,n=50),paste0(output,".significant_genes.txt"),col.names=NA)
    write.table(lrt$table,paste0(output,".all_genes.txt"),col.names=NA)
}

do.analysis1=function()
{

for (time_point in c(1,6,12))
{
    for (treatment1 in 1:4)
    {
        for (treatment2 in 1:4)
        {
            v_colnames=c()
            
            if (treatment1 < treatment2)
            {
                group=factor(c(1,1,1,2,2,2))
                for (r1 in 1:3)
                {
                  colname=paste0("X",time_point,"T",treatment1,"R",r1)
                  #print(colname)
                  v_colnames=append(v_colnames,colname)
                }
                # 1T4T1 is absent
                if (time_point !=1 | treatment2 != 4)
                {
                    colname=paste0("X",time_point,"T",treatment2,"R","1")
                    v_colnames=append(v_colnames,colname)
                }
                else
                {
                    group=factor(c(1,1,1,2,2))  
                }
                
                for (r1 in 2:3)
                {
                    colname=paste0("X",time_point,"T",treatment2,"R",r1)
                    #print(colname)
                    v_colnames=append(v_colnames,colname)
                }
                #print(v_colnames)
                x=all_counts[v_colnames]
                filename=paste0(time_point,".T",treatment1,"_vs_T",treatment2)
                de_comparison(x,group,filename)
            }
        }
    }
}
}

#time comparisons
do.analysis2=function()
{
  for (treatment in 1:4)
  {
    for (time_point1 in c(1,6,12))  
    {
      for (time_point2 in c(1,6,12))
      {
        v_colnames=c()
        
        if (time_point1 < time_point2)
        {
          group=factor(c(1,1,1,2,2,2))
          # 1T4T1 is absent
          if (time_point1 !=1 | treatment != 4)
          {
              colname=paste0("X",time_point1,"T",treatment,"R","1")
              v_colnames=append(v_colnames,colname)
          }
          else
          {
              group=factor(c(1,1,1,2,2))  
          }

          for (r1 in 2:3)
          {
            colname=paste0("X",time_point1,"T",treatment,"R",r1)
            #print(colname)
            v_colnames=append(v_colnames,colname)
          }          
          
          for (r1 in 1:3)
          {
              colname=paste0("X",time_point2,"T",treatment,"R",r1)
              #print(colname)
              v_colnames=append(v_colnames,colname)
          }
          
          print(v_colnames)
          x=all_counts[v_colnames]
          filename=paste0("X",time_point1,"_vs_X",time_point2,".T",treatment)
          de_comparison(x,group,filename)
        }
      }
    }
  }
}

do.analysis1()
do.analysis2()
