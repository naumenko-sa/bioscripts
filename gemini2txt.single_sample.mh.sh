#!/bin/bash
#   exports gemini.db database to gemini.db.txt file
#   database schema: https://gemini.readthedocs.io/en/latest/content/database_schema.html#the-variants-table
#   when using v.chr = g.chr AND v.gene = g.gene it becomes very slow
#   by default bcbio writes PASS only variants to the database

# for mh project - it does not have hgvsc/hgvsp

#PBS -l walltime=1:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

if [ -z $file ]
then
    file=$1
fi

sample=`gemini query -q "select name from samples" $file`

gemini query --header -q "select 
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
        v.aa_change as AA_change,
        v.codon_change as Codon_change,
	gts."$sample"
        from variants v, gene_detailed g
        where v.transcript=g.transcript and v.gene=g.gene and v.impact_severity <> 'LOW' and v.max_aaf_all < 0.01 and
          (v.gene='CACNA1S' or v.gene='RYR1' or v.gene='STAC3' or v.gene='TRDN' or v.gene='ASPH' or v.gene='JPH2' or 
          v.gene='CASQ1' or v.gene='ATP2A1' or v.gene='ATP2A2' or v.gene='CALM1' or v.gene='FKBP1A')" $file > ${file}.txt;
