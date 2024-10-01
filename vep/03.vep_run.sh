#!/bin/bash

sudo docker run -t -i -v /data/vep_data:/opt/vep/.vep:Z ensemblorg/ensembl-vep vep \
-o /opt/vep/.vep/variants.annotated.tsv \
--tab --cache  --merge --force_overwrite --offline --no_stats \
-input_file /opt/vep/.vep/variants.vep_input.txt \
--fasta /opt/vep/.vep/homo_sapiens/112_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
--symbol --numbers --biotype --total_length --canonical --hgvs --shift_hgvs 1 --af_gnomad \
--fields "Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,IMPACT,DISTANCE,STRAND,FLAGSSYMBOL,SYMBOL_SOURCE,HGNC_ID,BIOTYPE,CANONICAL,REFSEQ_MATCH,SOURCE,REFSEQ_OFFSET,GIVEN_REF,USED_REF,BAM_EDIT,EXON,INTRON,HGVSc,HGVSp,HGVS_OFFSET,gnomADe_AF,CLIN_SIG"
