#!/bin/bash
#   exports variant_effects from gemini.db database to gemini.db.variant_effects.txt file
#   database schema: https://gemini.readthedocs.io/en/latest/content/database_schema.html#the-variants-table
#   by default bcbio writes PASS only variants to the database
#   this is a version for refseq - some lof fields are absent

#PBS -l walltime=1:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

if [ -z $file ]
then
    file=$1
fi

sQuery="select 
	i.variant_id,
	i.gene,
	i.transcript,
	i.is_exonic,
	i.is_coding,
	i.is_lof,
	i.exon,
	i.codon_change,
	i.aa_change,
	i.aa_length,
	i.biotype,
	i.impact,
	i.impact_so,
	i.impact_severity,
	i.polyphen_pred,
	i.polyphen_score,
	i.sift_pred,
	i.sift_score,
	i.vep_canonical,
	i.vep_ccds,
	i.vep_hgvsc,
	i.vep_hgvsp 
	from variants v, variant_impacts i 
	where v.impact_severity <> 'LOW' and v.max_aaf_all < 0.01 and v.variant_id=i.variant_id"

echo $sQuery

gemini query --header -q "$sQuery" $file > ${file}.impacts.txt