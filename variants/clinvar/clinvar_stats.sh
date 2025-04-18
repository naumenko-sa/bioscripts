#!/bin/bash

# summarise all CLNSIG values
gunzip -c clinvar.vcf.gz | grep -v "^#" | awk -F  'CLNSIG=' '{print $2}' | awk -F ';' '{print $1}' | sort | uniq -c | awk '{print $2"\t"$1}' | sort -k2,2nr  > clinvar.vcf.gz.stat.tsv

# select only Pathogenic - 303083
# note filtering only by CLNSIG field
cat clinvar.vcf.gz.stat.tsv | grep -E "(Pathogenic|Likely_pathogenic)" | grep -v Conflicting_classifications_of_pathogenicity > clinvar.vcf.gz.stat.pathogenic.tsv

# parse CLNREFSTAT for Likely_pathogenic - 303083
# note filtering only by CLNSIG field, since Pathogenic occurs in other fields
gunzip -c clinvar.vcf.gz | grep -v "^#" | grep -E "CLNSIG=(Pathogenic|Likely_pathogenic)" | grep -v "CLNSIG=Conflicting_classifications_of_pathogenicity" | \
    awk '{print $8}' | awk -F 'CLNREVSTAT=' '{print $2}' | awk -F ';' '{print $1}' | sort | uniq -c | awk '{print $2"\t"$1}' | sort -k2,2nr > clinvar.vcf.gz.revstat.tsv
