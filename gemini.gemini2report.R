#Report for CHEO project according to project_cheo/report/cheo_report_format.doc

#moves last column to the position number
move_column=function (variants,number)
{
  columns = ncol(variants)
  if (number > 2)
    variants=variants[c(1:(number-1),columns,number:(columns-1))] 
  else
    variants=variants[c(1,columns,number:(columns-1))] 
  return(variants)
}

add_placeholder=function(variants,column_name,placeholder)
{
   variants[,column_name]=with(variants,placeholder)
   return(variants)
}

get_variants_from_db = function (dbname)
{
    dbname="NA12878-1-ensemble.db"
    con = dbConnect(RSQLite::SQLite(),dbname=dbname)
    dbListTables(con)

    qrySample="select name from samples"
    sample = dbGetQuery(con,qrySample)[1,1]

qryReport=paste0("select 
        v.ref as Ref,
        v.alt as Alt,
        v.impact as Variation,
        v.depth as Depth,
        v.qual_depth as Qual_depth,
        v.gene as Gene,
        g.ensembl_gene_id as Ensembl_gene_id,
        v.clinvar_disease_name as Clinvar,
        v.transcript as Ensembl_transcript_id,
        v.aa_length as AA_position,
        v.exon as Exon,
        v.pfam_domain as Pfam_domain,
        v.rs_ids as rsIDs,
        v.aaf_1kg_all as Maf_1000g,
        v.aaf_exac_all as Exac_maf,
        v.max_aaf_all as Maf_all,
        v.exac_num_het as Exac_het,
        v.exac_num_hom_alt as Exac_hom_alt,
        v.sift_score as Sift_score,
        v.polyphen_score as Polyphen_score,
        v.cadd_scaled as Cadd_score,gts,
        v.chrom as Chrom,
        v.start+1  as Pos,
        v.aa_change as AA_change
        from variants v, gene_detailed g
        where v.transcript=g.transcript and v.gene=g.gene and v.chrom = \"chr20\"");

        variants = dbGetQuery(con,qryReport)
    return (variants)
}

close_database = function()
{
  dbClearResult(variants)
  dbDisconnect(con)
}

#in the meanwhile using cheo.gemini2txt.single_sample.sh query in bash to decode 
#blob fields
get_variants_from_file = function (filename)
{
    variants = read.delim(filename, stringsAsFactors=FALSE)
    return(variants)
}

# return Hom, Het, or - if == hom_reference
genotype2zygocity = function (genotype_str,ref)
{
      #genotype_str = "A|A|B"
      #genotype_str = "./." - call not possible
      #genotype_str = "TCA/."
      #genotype_str = "G"
      #genotype_str="A/A"
      #greedy
      genotype_str = gsub("|","/",genotype_str,fixed=T)
      genotype_str = gsub("./.","Insufficient_coverage",genotype_str,fixed=T)
      #genotype_str = gsub(".","NO_CALL",genotype_str,fixed=T)
      
      if(grepl("Insufficient_coverage",genotype_str)){
          result = genotype_str
      }else{
          ar = strsplit(genotype_str,"/",fixed=T)
          len = length(ar[[1]])
          if (len == 2)
          {
            if (ar[[1]][1] == ar[[1]][2]){
              if (ar[[1]][1] == ref)
                  result = "-"
              else
                  result = "Hom"
            }else
              result = "Het"
          }else{
            result = genotype_str
          }
      }
      return(result)
}

#output : family.ensemble.txt
create_report = function(family,samples)
{
    #test: 3 samples in a family
    #family="166"
    #samples=c("166_3_5","166_4_10","166_4_8")
    #suffix = "gatk-haplotype"
    #suffix = "ensemble"
  
    #test: 1 sample in a familty
    #family="NA12878-1"
    #samples=c("NA12878.1")
  
    file=paste0(family,"-ensemble.db.txt")
  
    variants = get_variants_from_file(file)

  #field1 - Position
  variants$Position=with(variants,paste(Chrom,Pos,sep=':'))
  columns = ncol(variants)
  variants=variants[c(columns,1:columns-1)]

  #field2 - UCSC link
  sUCSC1="=HYPERLINK(\"http://genome.ucsc.edu/cgi-bin/hgTracks?hgt.out3=10x&position="
  sUCSC2="\",\"UCSC_link\""
  variants$UCSC_Link=with(variants,paste(sUCSC1,Position,sUCSC2,sep=''))

  #fields 3,4: Ref,Alt

  # field5 - Zygocity
  # use new loader vcf2db.py - with flag  to load plain text
  # for genotype and depth - Noah
  # otherwise have to decode BLOB 
  # snappy decompression
  # https://github.com/arq5x/gemini/issues/700
  # https://github.com/lulyon/R-snappy
  #variants = cbind(variants,lapply(variants$gts.100940,genotype2zygocity))
    for(sample in samples)
    {
      #DEBUG: gene = IL20RA
      #sample=samples[1]
      zygocity_column_name = paste0("Zygocity.",sample)
      #t = lapply(variants[,paste0("gts.",sample),"Ref"],genotype2zygocity)
      #t = lapply(variants[,paste0("gts.",sample),"Ref"],genotype2zygocity)
      t=unlist(mapply(genotype2zygocity,variants[,paste0("gts.",sample)],variants[,"Ref"]))
      variants[,zygocity_column_name] = unlist(t)
    
      burden_column_name = paste0("Burden.",sample)
      t = subset(variants, get(zygocity_column_name) == 'Hom' | get(zygocity_column_name) == 'Het',select=c("Ensembl_gene_id",zygocity_column_name))
      df_burden = count(t,'Ensembl_gene_id')    
      colnames(df_burden)[2] = burden_column_name
      variants = merge(variants,df_burden,all.x=T)
      variants[,burden_column_name][is.na(variants[,burden_column_name])] = 0
    }

    #field 6 - Variation
    #add splicing and splicing extended annotation - just use RefSeq annotation and chr:pos
    #in the meawhile we have splice_region_variant
    #http://sequenceontology.org/browser/current_svn/term/SO:0001630
    #1-3 bases of exon, 3-8 bases of intron

    #field 7 - Info, we want all transcripts, not it is more more severily affected
    #gemini bug: https://github.com/arq5x/gemini/issues/798
    ensembl_refseq = read.delim(paste0(reference_tables_path,"/ensembl_refseq_no_duplicates.txt"), stringsAsFactors=F)
    variants = merge(variants,ensembl_refseq,all.x=T)
    variants$Info = with(variants,paste(Gene,Refseq_mrna,Codon_change,Protein_change,sep=':'))

    #fields 8,9 - Depth, Qual_depth

    #field 10 - Alt_depth - from v.gt_alt_depth
    #when multiple callers used, AD is not set
    for(sample in samples)
    {
      new_name = paste0("Alt_depths.",sample)
      setnames(variants, paste0("gt_alt_depths.",sample),new_name)
      #variants[,new_name] = with(variants,gsub("-1","Multiple_callers",get(new_name),fixed=T))
    }

    #fields 11,12 - Gene, ENS_ID

    #field13 - Gene_description - from biomart
    gene_descriptions = read.delim2(paste0(reference_tables_path,"/ensembl_w_description.txt"), stringsAsFactors=F)
    variants = merge(variants,gene_descriptions,by.x = "Ensembl_gene_id",by.y = "ensembl_gene_id",all.x=T)

    #field14 - Omim_gene_description - from omim text file
    omim = read.delim2(paste0(reference_tables_path,"/omim.forannotation2"), stringsAsFactors=FALSE)
    variants = merge(variants,omim,all.x=T)

    #field15 - Omim_inheritance 
    #from Kristin xls
    omim_inheritance = read.delim(paste0(reference_tables_path,"/omim_inheritance.txt"), stringsAsFactors=F)
    variants = merge(variants,omim_inheritance,all.x=T)

    #field16 - Orphanet
    orphanet = read.delim(paste0(reference_tables_path,"/orphanet.deduplicated.txt"), stringsAsFactors=F)
    variants = merge(variants,orphanet,all.x=T)
    # mind the samples order: it will influence the Trio
    #fields17-18 Clinvar, Ensembl Transcript ID

    #fields19-20-21-22 - protein change, aa_position, exon, pfam domain

    #fields 23-24, will add when all samples will be done
    variants = add_placeholder(variants,"Frequency_in_C4R","Frequency_in_C4R")
    variants = add_placeholder(variants,"Seen_in_C4R_samples","Seen_in_C4R_samples")

    #field 25, 26 - rsIDs, Maf_1000g

    #fields 27-29 EVS frequencies

    #fields 30-31: Exac_maf, Maf_all
    #fields 32-33: Exac scores
    exac_scores = read.delim(paste0(reference_tables_path,"/exac_scores.txt"), stringsAsFactors=F)
    variants = merge(variants,exac_scores,all.x=T)

    #field 34-35  Exac het, exac_hom_alt

    #field36 - Conserved in 29 mammals instead of phastcons
    #https://www.biostars.org/p/150152/

    #fields 37-39: sift,polyphen,cadd

    #field40
    variants = add_placeholder(variants,"Trio_coverage","")
    n_sample = 1
    prefix = ""
  
    #order gts column in the same way as 
    variants$gts=""
    for(sample in samples)
    {
        column = paste0("gt_depths.",sample)
        if (n_sample>1) prefix="/"
        variants$Trio_coverage = with(variants,paste0(Trio_coverage,prefix,get(column)))
    
        column = paste0("gts.",sample)
        if (n_sample>1) prefix=","
            variants$gts = with(variants,paste0(gts,prefix,get(column)))
    
        n_sample = n_sample+1
    }

    #substitute -1 to 0
    for (field in c("EVS_maf_aa","EVS_maf_ea","EVS_maf_all","Maf_1000g","Exac_maf","Maf_all","Exac_het","Exac_hom_alt","Trio_coverage"))
    {
        variants[,field] = with(variants,gsub("-1","0",get(field),fixed=T))  
    }

    for (field in c(paste0("Alt_depths.",samples)))
    {
        variants[,field] = with(variants,gsub("-1",NA,get(field),fixed=T))  
    }


    #fields 41-42 - imprinting
    imprinting = read.delim(paste0(reference_tables_path,"/imprinting.txt"), stringsAsFactors=FALSE)
    variants = merge(variants,imprinting,all.x=T)

    #field 43 - pseudoautosomal
    pseudoautosomal = read.delim(paste0(reference_tables_path,"/pseudoautosomal.txt"), stringsAsFactors=F)
    variants = merge(variants,pseudoautosomal,all.x=T)

    select_and_write(variants,samples,paste0(family,".ensemble"))
}

#final selection and order
select_and_write = function(variants,samples,prefix)
{
    variants = variants[c(c("Position","UCSC_Link","Ref","Alt"),paste0("Zygocity.",samples),c("Gene"),
                        paste0("Burden.",samples),c("gts","Variation","Info","Depth","Quality"),
                        paste0("Alt_depths.",samples),c("Trio_coverage","Ensembl_gene_id","Gene_description","Omim_gene_description","Omim_inheritance",
                                                        "Orphanet", "Clinvar","Ensembl_transcript_id","Protein_change","AA_position","Exon","Pfam_domain",
                                                        "Frequency_in_C4R","Seen_in_C4R_samples","rsIDs","Maf_1000g","EVS_maf_aa","EVS_maf_ea","EVS_maf_all",
                                                        "Exac_maf","Maf_all", "Exac_pLi_score","Exac_missense_score","Exac_het","Exac_hom_alt",
                                                        "Conserved_in_29_mammals","Sift_score","Polyphen_score","Cadd_score",
                                                        "Imprinting_status","Imprinting_expressed_allele","Pseudoautosomal"))]
  
    write.table(variants,paste0(prefix,".txt"),quote=F,sep = ";",row.names=F)  
}

# merges ensembl and gatk-haplotype reports to 
# - fill alt_depths, Trio_coverage, qual_depths columns
# fix NO_CALL issue
#   complex alleles like T > C,TA are divided into two variants : T>C, T>TA
#   why are they not separated before loading to gemini - check
#   indels called differently should be reported from GATK
merge_reports = function(family,samples)
{
    #setwd("/home/sergey/Desktop/project_cheo/2016-11-09_rerun10")
    #family = "166"
    # mind the samples order: it will influence the Trio
    #samples=c("166_3_5","166_4_10","166_4_8")
  
    ensemble_file = paste0(family,".ensemble.txt")
    
    ensemble = read.csv(ensemble_file, sep=";", quote="", stringsAsFactors=F)
    ensemble$superindex=with(ensemble,paste(Position,Ref,Alt,sep='-'))
    
    gatk_file = paste0(family,"-gatk-haplotype.decomposed.table")
    gatk = read.delim(gatk_file, stringsAsFactors=F)
    gatk$superindex=with(gatk,paste(paste0("chr",CHROM,":",POS),REF,ALT,sep='-'))
    gatk[c("CHROM","POS","REF","ALT")]=NULL
    
    ensemble = merge(ensemble,gatk,by.x = "superindex", by.y="superindex",all.x = T)
    
    ensemble$Depth = ensemble$DP
    n_sample = 1
    prefix = ""
    ensemble$Trio_coverage=""
    
    for(sample in samples)
    {
        #R fixes numerical column names with X
        column = paste0("X",sample,".DP")
        if (n_sample>1) prefix="/"
        ensemble$Trio_coverage = with(ensemble,paste0(Trio_coverage,prefix,get(column)))
      
        column = paste0("Alt_depths.",sample)
        column_gatk = paste0("X",sample,".AD")
        
        ensemble[,column] = ensemble[,column_gatk]
      
        n_sample = n_sample+1
    }
    
    for (i in 1:nrow(ensemble))
    {
        for (sample in samples)
        {
            field = paste0("Alt_depths.",sample)
            
            ensemble[i,field]=strsplit(ensemble[i,field],",",fixed=T)[[1]][2]
        }
    }
    
    ensemble[c("DP",paste0("X",samples,".DP"),paste0("X",samples,".AD"))]=NULL
    #ensemble[c("DP",paste0(samples,".DP"),paste0(samples,".AD"))]=NULL
    
    freebayes_file = paste0(family,"-freebayes.decomposed.table")
    freebayes = read.delim(freebayes_file, stringsAsFactors=F)
    freebayes$superindex=with(freebayes,paste(paste0("chr",CHROM,":",POS),REF,ALT,sep='-'))
    freebayes[c("CHROM","POS","REF","ALT")]=NULL
    ensemble = merge(ensemble,freebayes,by.x = "superindex", by.y="superindex",all.x = T)
    
    for (i in 1:nrow(ensemble))
    {
        #if (ensemble[i,"Trio_coverage"]=="NA/NA/NA")
        if(grepl("NA",ensemble[i,"Trio_coverage"]))
        {
            ensemble[i,"Depth"] = ensemble[i,"DP"]
            for (sample in samples)
            {
                field_depth = paste0("Alt_depths.",sample)
                field_bayes = paste0("X",sample,".AO")
                #field_bayes = paste0(sample,".AO")
                      
                ensemble[i,field_depth] = ensemble[i,field_bayes]
        
            }
            n_sample = 1
            prefix = ""
            ensemble[i,"Trio_coverage"]=""
            
            for(sample in samples)
            {
                column = paste0("X",sample,".DP")
                if (n_sample>1) prefix="/"
                ensemble[i,"Trio_coverage"] = paste(ensemble[i,"Trio_coverage"],ensemble[i,column],sep = prefix)
              
                n_sample = n_sample+1
            }
        }
    }
    
    ensemble[c("DP",paste0("X",samples,".DP"),paste0("X",samples,".AO"))]=NULL
    
    platypus_file = paste0(family,"-platypus.decomposed.table")
    platypus = read.delim(platypus_file, stringsAsFactors=F)
    platypus$superindex=with(platypus,paste(paste0("chr",CHROM,":",POS),REF,ALT,sep='-'))
    platypus[c("CHROM","POS","REF","ALT")]=NULL
    ensemble = merge(ensemble,platypus,by.x = "superindex", by.y="superindex",all.x = T)
    
    for (i in 1:nrow(ensemble))
    {
      if(grepl("NA",ensemble[i,"Trio_coverage"]))
      #if (ensemble[i,"Trio_coverage"]=="NA/NA/NA")
      {
        ensemble[i,"Depth"] = ensemble[i,"TC"]
        for (sample in samples)
        {
          field_depth = paste0("Alt_depths.",sample)
          field_bayes = paste0("X",sample,".NV")
          
          #sometimes freebayes has 10,10,10 for decomposed alleles
          ensemble[i,field_depth] = strsplit(ensemble[i,field_bayes],",",fixed=T)[[1]][1]
          
        }
        n_sample = 1
        prefix = ""
        ensemble[i,"Trio_coverage"]=""
        
        for(sample in samples)
        {
          column = paste0("X",sample,".NR")
          if (n_sample>1) prefix="/"
          #sometimes freebayes has 10,10,10 for decomposed alleles
          cov_value=strsplit(ensemble[i,column],",",fixed=T)[[1]][1]
          ensemble[i,"Trio_coverage"] = paste(ensemble[i,"Trio_coverage"],cov_value,sep = prefix)
          
          n_sample = n_sample+1
        }
      }
    }
    
    
    ensemble[c("TC",paste0("X",samples,".NV"),paste0("X",samples,".NR"))]=NULL
    #ensemble[c("TC",paste0(samples,".NV"),paste0(samples,".NR"))]=NULL
    ensemble[,"Trio_coverage"] = with(ensemble,gsub("NA","0",get("Trio_coverage"),fixed=T))  
   
    for (i in 1:nrow(ensemble))
    {
        if (is.na(ensemble[i,"Depth"]))
        {
            l=strsplit(ensemble[i,"Trio_coverage"],"/")[[1]]
            ensemble[i,"Depth"]=sum(as.integer(l))
        }
        for (sample in samples)
        {
            field_depth = paste0("Alt_depths.",sample)
            if (is.na(ensemble[i,field_depth]))
                ensemble[i,field_depth]=0
        }
    }
    
    select_and_write(ensemble,samples,family)
}

library(RSQLite)
library(stringr)
library(genetics)
library(data.table)
library(plyr)

reference_tables_path="~/Desktop/reference_tables"

#Hernan samples
setwd("/home/sergey/Desktop/project_muscular/Muscle2/")
family="muscle2"
samples=c("Muscle2_filtered")
create_report(family,samples)

setwd("/home/sergey/Desktop/project_muscular/DMD/")
family="dmd"
samples=c("DMD")
create_report(family,samples)

#cheo 10 samples
setwd("/home/sergey/Desktop/project_cheo/2016-11-09_rerun10")
family="166"
samples=c("166_3_5","166_4_10","166_4_8")
create_report(family,samples)
merge_reports(family,samples)

family="181"
samples = c("181_121141J","181_WG0927")
create_report(family,samples)
merge_reports(family,samples)

family="241"
samples=c("241_44845","241_52062","241_52063")
create_report(family,samples)
merge_reports(family,samples)

family = "246"
samples = c("246_90137","246_CH0015","246_CH0016")
create_report(family,samples)
merge_reports(family,samples)

family="380"
samples=c("380_120890B","380_120891B")
create_report(family,samples)
merge_reports(family,samples)

family="391"
samples=c("391_121030T","391_121031C","391_CH0073")
create_report(family,samples)
merge_reports(family,samples)

family="394"
samples=c("394_60638BD")
create_report(family,samples)
merge_reports(family,samples)

family="411"
samples=c("411_G0071AG","411_G0091AG")
create_report(family,samples)
merge_reports(family,samples)

family="412"
samples=c("412_120880N","412_120886B","412_120887D")
create_report(family,samples)
merge_reports(family,samples)

family="417"
samples=c("417_120882D")
create_report(family,samples)
merge_reports(family,samples)


#mutations for katie
setwd("/home/sergey/Desktop/project_katie_csc_large")
families=c("wG432_DMSO","wG432_LGK","wG440_DMSO","wG440_LGK","wG472_DMSO","wG472_LGK",
           "wG481_DMSO","wG481_LGK","wG510_DMSO","wG510_LGK","wG511_DMSO","wG511_LGK",
           "wG523_DMSO","wG523_LGK","wG564_DMSO","wG564_LGK")
for (family in families)
{
  samples=c("gbm")
  create_report(family,samples)
}

#V
setwd("/home/sergey/Desktop/project_exomes/1_v/2016-12-19_new_report_decomposed/")
family="bnu7823"
samples=c("nu7823")
create_report(family,samples)
merge_reports(family,samples)

#N
setwd("/home/sergey/Desktop/project_exomes/2_n/2016-12-19_new_report_decomposed")
family="b100940"
samples=c("100940")
create_report(family,samples)
merge_reports(family,samples)

#mh
setwd("/home/sergey/Desktop/project_mh/B175/")
family="B175"
samples = c("1130-BD-B175","2064-BA-B175")
create_report(family,samples)
merge_reports(family,samples)

# R substitutes - with . in sample names in columns
setwd("/home/sergey/Desktop/project_cheo/2016-12-26_reports_50_families")
families <- unlist(read.table("families_ready.txt", quote="\"", comment.char="", stringsAsFactors=FALSE))

for (family in families)
{
    setwd(family)
    samples = unlist(read.table("samples.txt", quote="\"", comment.char="", stringsAsFactors=FALSE))
    samples = gsub("-",".",samples)
    create_report(family,samples)
    merge_reports(family,samples)
    setwd("..")
}
