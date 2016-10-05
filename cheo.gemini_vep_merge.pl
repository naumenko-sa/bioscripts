#!/usr/bin/perl

#merges results from gemini2txt and vcfvep2txt in 1 file 
#quick dirty solution should be possible to extract everything from gemini

#vcf FIELDs CHR\tPOS\tID\tREF\tALT\tGT\tConsequence\tCodons\tAmino_acids\tGene\tSYMBOL\tFeature\tEXON\tPolyPhen\tSIFT\tProtein_position\tBIOTYPE\tCANONICAL\tCCDS\tLoF\tLoF_filter\tLoF_flags"
#vcf KEEP: GT-no6 canonical-18 ccds-19

#gemini FIELDs
#chrom	start0based	start1based	end1based	ref	alt	qual	type	sub_type	gene	
#ensembl_gene_id	transcript	is_exonic	is_coding	is_lof is_splicing	exon	codon_change	aa_change	biotype	
#impact	impact_so	impact_severity	polyphen_pred	polyphen_score sift_pred	sift_score	pfam_domain	depth	in_hom_run	qual_depth	
#in_dbsnp	rs_ids	in_1kg	in_exac	aaf_exac_all max_aaf_all	exac_num_het	exac_num_hom_alt	exac_num_chroms	in_omim	clinvar_causal_allele	clinvar_sig	clinvar_disease_name	
#exome_chip	cyto_band	rmsk	in_cpg_island	in_segdup	is_conserved	recomb_rate	cadd_raw	cadd_scaled	grc	gms_illumina	in_cse	info	
#exac_pLi	exac_mis_z	old_call	omim

open(VEP,$ARGV[0]);
open(GEMINI,$ARGV[1]);

%vep_db;
%gemini_db;

while(<GEMINI>)
{
    chomp;
    @ar=split("\t",$_);
    $ar[0] =~ s/chr//;
    $id=$ar[0]."-".$ar[2]."-".$ar[4]."-".$ar[5];
    $gemini_db{$id}=$ar[35]."\t".$ar[36];
}

while(<VEP>)
{
    chomp;
    @ar=split("\t",$_);
    $id=$ar[0]."-".$ar[1]."-".$ar[3]."-".$ar[4];
    if ($ar[0] == "CHR")
    {
	print $_."\tAAF_EXAC_ALL\tMAX_AAF_ALL\n";
    }
    else
    {
	print $_."\t".$gemini_db{$id}."\n";
    }
}


close(VEP);
close(GEMINI);