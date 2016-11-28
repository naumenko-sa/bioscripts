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

add_placeholder=function(variants,column_name,placeholder,number)
{
   variants[,column_name]=with(variants,placeholder)
   return(move_column(variants,number))
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

# return Homozygous Heterozygous or Multiple het
genotype2zygocity = function (genotype_str)
{
      genotype_str = "A|A|B"
      #genotype_str = "./."
      #genotype_str = "TCA/."
      #genotype_str = "G"
      #greedy
      genotype_str = gsub("|","/",genotype_str,fixed=T)
      
      ar = strsplit(genotype_str,"/",fixed=T)
      result = genotype_str
      return(result)
}

create_report = function(family,samples)
{
  #file="417-ensemble.db.txt"
  #sample="417_120882D"
  
  samples=c("166_3_5","166_4_10","166_4_8")
  family="166"
  
  samples=c("100940")
  family = "b100940"
  
  file=paste0(family,"-ensemble.db.txt")
  
  variants = get_variants_from_file(file)


#field1 - Position
variants$Position=with(variants,paste(Chrom,Pos,sep=':'))
columns = ncol(variants)
variants=variants[c(columns,1:columns-1)]

#field2 - UCSC link
sUCSC1="=HYPERLINK(\"http://genome.ucsc.edu/cgi-bin/hgTracks?hgt.out3=10x&position="
sUCSC2="\",\"UCSC link\""
variants$UCSC_Link=with(variants,paste(sUCSC1,Position,sUCSC2,sep=''))
variants=move_column(variants,2)

#fields 3,4

# field5 - Zygocity - 
# use new loader vcf2db.py - with flag  to load plain text
# for genotype and depth - Noah
# otherwise have to decode BLOB 
# snappy decompression
# https://github.com/arq5x/gemini/issues/700
# https://github.com/lulyon/R-snappy
variants=add_placeholder(variants,new_name,"test",7)
#variants = cbind(variants,lapply(variants$gts.100940,genotype2zygocity))
for(sample in samples)
{
    new_name = paste0("Zygocity.",sample)
    #setnames(variants, paste0("gts.",sample),new_name)
    variants[,new_name] = with(variants,genotype2zygocity(paste0("gts.",sample)))
}

#field 6 - Variation
#add splicing and splicing extended annotation - just use RefSeq annotation and chr:pos
#in the meawhile we have splice_region_variant
#http://sequenceontology.org/browser/current_svn/term/SO:0001630
#1-3 bases of exon, 3-8 bases of intron

#field 7 - Info
#variants=add_placeholder(variants,"Info",paste("Gene","NCBI_TRANSCRIPT","NUC_CHANGE","PROT_CHANGE",sep=':'),7)
ensembl_refseq = read.delim(paste0(reference_tables_path,"/ensembl_refseq.txt"), stringsAsFactors=FALSE)
#variants = merge(variants,ensembl_refseq,all.x=T)


#variants = cbind(variants,str_split_fixed(variants$Protein_change,":",2))
#variants = rename(variants,c("1"="Junk1","Protein_change"="Junk2","2"="Protein_change"))

#variants = cbind(variants,str_split_fixed(variants$Codon_change,":",2))
#variants = rename(variants,c("1"="Junk3","2"="Info_codon_change"))

#variants$Info = with(variants,paste(Gene,Refseq_mrna,Info_codon_change,Protein_change,sep=':'))

#fields 8,9 - Depth, Qual_depth

#field 10 - Alt_depth - from v.gt_alt_depth
#alt_column_name = paste0("gt_alt_depths.",sample)
#vriants = rename(variants,replace=c(column_name="Alt_depth"))

#fields 11,12 - Gene, ENS_ID

#field13 - from biomart
#variants=add_placeholder(variants,"Gene_description","Gene_description",13)
gene_descriptions = read.delim2(paste0(reference_tables_path,"/ensembl_w_description.txt"), stringsAsFactors=FALSE)
variants = merge(variants,gene_descriptions,by.x = "Ensembl_gene_id",by.y = "ensembl_gene_id",all.x=T)

#field14 - from omim text file
#variants = add_placeholder(variants,"Omim_gene_description","Omim_gene_description",14)
omim = read.delim2(paste0(reference_tables_path,"/omim.forannotation2"), stringsAsFactors=FALSE)
variants = merge(variants,omim,all.x=T)

#field15 - from Kristin xls
#variants = add_placeholder(variants,"Omim_inheritance","Omim_inheritance",15)
omim_inheritance = read.delim(paste0(reference_tables_path,"/omim_inheritance.txt"), stringsAsFactors=F)
variants = merge(variants,omim_inheritance,all.x=T)

#field16 - Orphanet
variants = add_placeholder(variants,"Orphanet","Orphanet",16)

#fields17-18

#field19 - protein change
#variants = add_placeholder(variants,"Protein_change","Protein_change",19)
#variants$Protein_change = with(variants,paste("p.",AA_change,AA_position,sep=' '))

#fields 23-24
variants = add_placeholder(variants,"Frequency_in_C4R","Frequency_in_C4R",23)
variants = add_placeholder(variants,"Seen_in_C4R_samples","Seen_in_C4R_samples",24)

#fields 27-28
variants = add_placeholder(variants,"EVS_maf","EVS_maf",27)
variants = add_placeholder(variants,"EVS_genotype_counts","EVS_genotype_counts",28)

#fields 31-32
exac_scores = read.delim(paste0(reference_tables_path,"/exac_scores.txt"), stringsAsFactors=F)
#variants = add_placeholder(variants,"Exac_pLi_score","Exac_pLi_score",31)
#variants = add_placeholder(variants,"Exac_missense_score","Exac_missense_score",32)
variants = merge(variants,exac_scores,all.x=T)

#field 35
#https://www.biostars.org/p/150152/
variants = add_placeholder(variants,"Phast_cons_score","Phast_cons_score",35)

#fields39-42
variants = add_placeholder(variants,"Trio_coverage","Trio_coverage",39)

#variants = add_placeholder(variants,"Imprinting_status","Imprinting_status",40)
#variants = add_placeholder(variants,"Imprinting_expressed_allele","Imprinting_expressed_allele",41)
imprinting = read.delim(paste0(reference_tables_path,"/imprinting.txt"), stringsAsFactors=FALSE)
variants = merge(variants,imprinting,all.x=T)

variants = add_placeholder(variants,"Pseudoautosomal","Pseudoautosomal",42)

#full set
#variants = variants[c(c("Position","UCSC_Link","Ref","Alt","Zygocity","Variation","Info","Depth","Qual_depth"),paste0("gts.",samples),
#                      c("Gene","Ensembl_gene_id","Gene_description","Omim_gene_description","Omim_inheritance","Orphanet",
#                      "Clinvar","Ensembl_transcript_id","Protein_change","AA_position","Exon","Pfam_domain","Frequency_in_C4R",
#                      "Seen_in_C4R_samples","rsIDs","Maf_1000g","EVS_maf","EVS_genotype_counts","Exac_maf","Maf_all",
#                      "Exac_pLi_score","Exac_missense_score","Exac_het","Exac_hom_alt","Phast_cons_score","Sift_score",
#                      "Polyphen_score","Cadd_score","Trio_coverage","Imprinting_status","Imprinting_expressed_allele",
#                      "Pseudoautosomal"))]

#do not report placeholders
variants = variants[c(c("Position","UCSC_Link","Ref","Alt","Zygocity","Variation","Depth"),paste0("gts.",samples),
                      c("Gene","Ensembl_gene_id","Gene_description","Omim_gene_description","Omim_inheritance",
                        "Clinvar","Ensembl_transcript_id","Protein_change","AA_position","Exon","Pfam_domain",
                        "rsIDs","Maf_1000g","Exac_maf","Maf_all",
                        "Exac_pLi_score","Exac_missense_score","Exac_het","Exac_hom_alt","Sift_score",
                        "Polyphen_score","Cadd_score","Imprinting_status","Imprinting_expressed_allele"))]


write.table(variants,paste0(family,".txt"),quote=F,sep = ";",row.names=F)  
#close_database()
}

setwd("~/Desktop/project_cheo/2016-10-28_gemini_test/")

reference_tables_path="~/Desktop/reference_tables"
library(RSQLite)
library(stringr)
library(genetics)
library(data.table)

create_report("166",c("166_3_5","166_4_10","166_4_8"))
create_report("181",c("181_121141J","181_WG0927"))
create_report("241",c("241_44845","241_52062","241_52063"))
create_report("246",c("246_90137","246_CH0015","246_CH0016"))
create_report("380",c("380_120890B","380_120891B"))
create_report("391",c("391_121031C","391_CH0073"))
create_report("394",c("394_60638BD"))
create_report("411",c("411_G0071AG","411_G0091AG"))
create_report("412",c("412_120880N","412_120886B","412_120887D"))
create_report("417",c("417_120882D"))

setwd("~/Desktop/project_exomes/2_n/2016-11-17_ensemble_calls/")
setwd("~/Desktop/project_exomes")
create_report("b100940",c("100940"))