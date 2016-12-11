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

test = function()
{
  library(stringr)
  samples=c("166_3_5","166_4_10","166_4_8")
  family="166"
  file=paste0(family,"-ensemble.db.txt")
  variants = get_variants_from_file(file)
  variants$allele_pool = paste(variants$gts.166_4_8,variants$gts.166_4_10)
  variants$allele_pool = gsub("/"," ",variants$test,fixed=T)
  variants$allele_pool = gsub("|"," ",variants$test,fixed=T)
  variants$allele_pool = gsub("."," ",variants$test,fixed=T)
  #variants$allele_pool = gsub("  "," ",variants$test,fixed=T)
  
  variants$allele1 = str_split_fixed(variants$gts.166_3_5,"/",2)
  
  variants$inherited = with(variants,{t=strsplit(allele_pool," ",fixed=T);t[2]})
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

#suffix = [ensemble | gatk-haplotype]
create_report = function(family,samples,suffix)
{
  #test: 3 samples in a family
  #family="166"
  samples=c("166_3_5","166_4_10","166_4_8")
  #suffix = "gatk-haplotype"
  #suffix = "ensemble"
  
  #test: 1 sample in a familty
  #family="NA12878-1"
  #samples=c("NA12878.1")
  
  file=paste0(family,"-",suffix,".db.txt")
  
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
for(sample in samples)
{
  column = paste0("gt_depths.",sample)
  if (n_sample>1) prefix="/"
  variants$Trio_coverage = with(variants,paste0(Trio_coverage,prefix,get(column)))
  n_sample = n_sample+1
}
#substitute -1 to 0
for (field in c("EVS_maf_aa","EVS_maf_ea","EVS_maf_all","Maf_1000g","Exac_maf","Maf_all","Exac_het","Exac_hom_alt","Trio_coverage"))
{
  variants[,field] = with(variants,gsub("-1","0",get(field),fixed=T))  
}

for (field in c(paste0("Alt_depths.",samples)))
{
  variants[,field] = with(variants,gsub("-1","0",get(field),fixed=T))  
}


#fields 41-42 - imprinting
imprinting = read.delim(paste0(reference_tables_path,"/imprinting.txt"), stringsAsFactors=FALSE)
variants = merge(variants,imprinting,all.x=T)

#field 43 - pseudoautosomal
pseudoautosomal = read.delim(paste0(reference_tables_path,"/pseudoautosomal.txt"), stringsAsFactors=F)
variants = merge(variants,pseudoautosomal,all.x=T)

    select_and_write(variants,samples,paste0(family,".",suffix))
}

#final selection and order
select_and_write = function(variants,samples,prefix)
{
    variants = variants[c(c("Position","UCSC_Link","Ref","Alt"),paste0("Zygocity.",samples),c("Gene"),
                        paste0("Burden.",samples),c("gts","Variation","Info","Depth","Qual_depth"),
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
    setwd("/home/sergey/Desktop/project_cheo/2016-11-09_rerun10")
    family = "166"
    # mind the samples order: it will influence the Trio
    samples=c("166_3_5","166_4_10","166_4_8")
    ensemble_file = paste0(family,".ensemble.txt")
    gatk_file = paste0(family,".gatk-haplotype.txt")
    
    ensemble = read.csv(ensemble_file, sep=";", quote="", stringsAsFactors=F)
    gatk = read.csv(gatk_file, sep=";", quote="", stringsAsFactors=F)
    gatk.depths = gatk[c("Position","Ref","Alt",paste0("Alt_depths.",samples),"Trio_coverage","Qual_depth")]

    ensemble[c(paste0("Alt_depths.",samples),"Trio_coverage","Qual_depth")]=NULL
    
    ensemble$superindex=with(ensemble,paste(Position,Ref,Alt,sep='-'))
    gatk.depths$superindex=with(gatk.depths,paste(Position,Ref,Alt,sep='-'))
    gatk.depths[c("Ref","Alt","Position")]=NULL
    
    ensemble = merge(ensemble,gatk.depths,by.x = "superindex", by.y="superindex",all.x = T)
    
    freebayes = read.delim("166-freebayes.decomposed.table", stringsAsFactors=F)
    freebayes$superindex=with(freebayes,paste(paste0("chr",CHROM,":",POS),REF,ALT,sep='-'))
    freebayes[c("CHROM","POS","REF","ALT")]=NULL
    ensemble = merge(ensemble,freebayes,by.x = "superindex", by.y="superindex",all.x = T)
    
    platypus = read.delim("166-platypus.decomposed.table", stringsAsFactors=F)
    platypus$superindex=with(platypus,paste(paste0("chr",CHROM,":",POS),REF,ALT,sep='-'))
    ensemble = merge(ensemble,platypus,by.x = "superindex", by.y="superindex",all.x = T)
    
    #if alt_depth and trio_depth is absent, get from freebayes or platypus
    
    for (i in 1:nrow(ensemble))
    {
        for (sample in samples)
        {
            field_depth = paste0("Alt_depths.",sample)
            field_bayes = paste0("X",sample,".AO")
            field_platypus = paste0("X",sample,".NV")
        
            if (is.na(ensemble[i,field_depth])) 
              ensemble[i,field_depth] = ensemble[i,field_bayes]
        
            if (is.na(ensemble[i,field_depth])) 
              ensemble[i,field_depth] = ensemble[i,field_platypus]
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
#setwd("/home/sergey/Desktop/project_cheo/2016-10-28_gemini_test")
setwd("/home/sergey/Desktop/project_cheo/2016-11-09_rerun10/ensemble_txt")
setwd("/home/sergey/Desktop/project_cheo/2016-11-09_rerun10/haplotype_txt")

create_report("NA12878-1",samples=c("NA12878.1"))

suffix="ensemble"
suffix="gatk-haplotype"

create_report("166",c("166_3_5","166_4_10","166_4_8"),suffix)
create_report("181",c("181_121141J","181_WG0927"),suffix)
create_report("241",c("241_44845","241_52062","241_52063"),suffix)
create_report("246",c("246_90137","246_CH0015","246_CH0016"),suffix)
create_report("380",c("380_120890B","380_120891B"),suffix)
create_report("391",c("391_121030T","391_121031C","391_CH0073"),suffix)
create_report("394",c("394_60638BD"),suffix)
create_report("411",c("411_G0071AG","411_G0091AG"),suffix)
create_report("412",c("412_120880N","412_120886B","412_120887D"),suffix)
create_report("417",c("417_120882D"),suffix)



