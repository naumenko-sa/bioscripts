#!/bin/bash

#exports gemini database to txt file to import to xls
#database schema: https://gemini.readthedocs.io/en/latest/content/database_schema.html#the-variants-table
#reports PASS only variants to the database

bname=`echo $1 | sed s/.db//`;
sample=$2

gemini query --header -q "select v.chrom,
			  v.start as start0based,
			  v.start+1  as start1based,
			  v.end as end1based,
			  v.ref,
			  v.alt,
			  v.qual,
			  v.type,
			  v.sub_type,
			  (v.gts)."$sample" as genotype,
			  v.gene,
			  g.ensembl_gene_id,
			  v.transcript,
			  v.is_exonic,
			  v.is_coding,
			  v.is_lof,
			  v.is_splicing,
			  v.exon,
			  v.codon_change,
			  v.aa_change,
			  v.biotype,
			  v.impact,
			  v.impact_so,
			  v.impact_severity,
			  v.polyphen_pred,
			  v.polyphen_score,
			  v.sift_pred,
			  v.sift_score,
			  v.pfam_domain,
			  v.depth,
			  v.in_hom_run,
			  v.qual_depth,
			  v.in_dbsnp,
			  v.rs_ids,
			  v.in_1kg,
			  v.in_exac,
			  v.aaf_exac_all,
			  v.max_aaf_all,
			  v.exac_num_het,
			  v.exac_num_hom_alt,
			  v.exac_num_chroms,
			  v.in_omim,
			  v.clinvar_causal_allele,
			  v.clinvar_sig,
			  v.clinvar_disease_name,
			  v.exome_chip,
			  v.cyto_band,
			  v.rmsk,
			  v.in_cpg_island,
			  v.in_segdup,
			  v.is_conserved,
			  v.recomb_rate,
			  v.cadd_raw,
			  v.cadd_scaled,
			  v.grc,
			  v.gms_illumina,
			  v.in_cse,
			  v.info
			  from variants v, gene_summary g
			  where v.chrom=g.chrom AND v.gene=g.gene" $1 > $bname.tmp;

head -n1 $bname.tmp | awk '{for (i=1;i<=NF;i++) printf $i"\t";printf "exac_pLi\texac_mis_z\told_call\tomim\n"}' > $bname.txt
cat $bname.tmp | sed 1d > $bname.body;
cheo.add_exac_scores.pl $bname.body > $bname.body1;
cheo.add_omim_string.pl $bname.body1 >>  $bname.txt
#rm $bname.tmp $bname.body $bname.body1