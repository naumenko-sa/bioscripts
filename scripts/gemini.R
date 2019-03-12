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