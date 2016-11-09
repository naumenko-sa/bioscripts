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

#in the meanwhile using gemini query in bash to decode 
#blob fields
get_variants_from_file = function (filename)
{
    variants = read.delim(filename, stringsAsFactors=FALSE)
    return(variants)
}

setwd("~/Desktop/project_cheo/2016-10-28_gemini_test/")
library(RSQLite)
library(stringr)

#variants = get_variants_from_db("NA12878-1-ensemble.db.txt")
variants = get_variants_from_file("NA12878-1-ensemble.db.txt")

#field1 - Position
variants$Position=with(variants,paste(Chrom,Pos,sep=':'))
columns = ncol(variants)
variants=variants[c(columns,1:columns-1)]

#field2 - UCSC link
sUCSC1="=HYPERLINK(\"http://genome.ucsc.edu/cgi-bin/hgTracks?hgt.out3=10x&position=\""
sUCSC2=",\"UCSC link\""
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
variants=add_placeholder(variants,"Zygocity","Zygocity",5)
variants$Zygocity = with(variants,gts)

#field 6 - Variation
#add splicing and splicing extended annotation - just use RefSeq annotation and chr:pos
#in the meawhile we have splice_region_variant
#http://sequenceontology.org/browser/current_svn/term/SO:0001630
#1-3 bases of exon, 3-8 bases of intron

#field 7 - Info
variants=add_placeholder(variants,"Info",paste("Gene","NCBI_TRANSCRIPT","NUC_CHANGE","PROT_CHANGE",sep=':'),7)
variants$Info = with(variants,paste(Gene,"NCBI_TRANSCRIPT","NUC_CHANGE","PROT_CHANGE",sep=':'))

#fields 8,9 - Depth, Qual_depth

#field 10 - Alt_depth - from v.gt_alt_depth
variants=add_placeholder(variants,"Alt_depth","Alt_depth",10)

#fields 11,12 - Gene, ENS_ID

#field13 - from biomart
variants=add_placeholder(variants,"Gene_description","Gene_description",13)
gene_descriptions = read.delim2("ensembl_w_description.txt", stringsAsFactors=FALSE)
variants = merge(variants,gene_descriptions,by.x = "Ensembl_gene_id",by.y = "ensembl_gene_id",all.x=T)
variants$Gene_description = with(variants,description)
variants = within(variants,rm(description))

#field14 - from omim text file
#variants = add_placeholder(variants,"Omim_gene_description","Omim_gene_description",14)
omim = read.delim2("omim.forannotation1", stringsAsFactors=FALSE)
variants = merge(variants,omim,all.x=T)

#field15 - from Kristin xls
variants = add_placeholder(variants,"Omim_inheritance","Omim_inheritance",15)
omim_inheritance = read.delim("omim_inheritance.txt", stringsAsFactors=FALSE)
variants = merge(variants,omim_inheritance,all.x=T)

#field16 - Orphanet
variants = add_placeholder(variants,"Orphanet","Orphanet",16)

#fields17-18

#field19 - protein change
variants = add_placeholder(variants,"Protein_change","Protein_change",19)
variants$Protein_change = with(variants,paste("p.",AA_change,AA_position,sep=' '))

#fields 23-24
variants = add_placeholder(variants,"Frequency_in_C4R","Frequency_in_C4R",23)
variants = add_placeholder(variants,"Seen_in_C4R_samples","Seen_in_C4R_samples",24)

#fields 27-28
variants = add_placeholder(variants,"EVS_maf","EVS_maf",27)
variants = add_placeholder(variants,"EVS_genotype_counts","EVS_genotype_counts",28)

#fields 31-32
variants = add_placeholder(variants,"Exac_pLi_score","Exac_pLi_score",31)
variants = add_placeholder(variants,"Exac_missense_score","Exac_missense_score",32)

#field 35
variants = add_placeholder(variants,"Phast_cons_score","Phast_cons_score",35)

#fields39-42
variants = add_placeholder(variants,"Trio_coverage","Trio_coverage",39)
variants = add_placeholder(variants,"Imprinting_status","Imprinting_status",40)
variants = add_placeholder(variants,"Imprinting_expressed_allele","Imprinting_expressed_allele",41)
variants = add_placeholder(variants,"Pseudoautosomal","Pseudoautosomal",42)

#close_database()
