#!/bin/bash

# run orthomcl to build cogs for a group of species in species
# manual: http://orthomcl.org/common/downloads/software/v2.0/UserGuide.txt

#1.copy files - use 1 isoform per locus, cd-hit
input_dir="/mnt/lustre/snaumenko1/project_gammarus/3cds/cdhit_n_longest_cds_per_locus"

suffix="cd-hit.fixed"

for f in `cat species`;
do
    cp $input_dir/${f}.$suffix.cds.fasta .
    cp $input_dir/${f}.$suffix.pep.fasta .
done

#2.remove . from identifiers
cat gam8.4.cd-hit.fixed.cds.fasta | sed s/"gam8.4"/"gam8_4"/ > gam8_4.cd-hit.fixed.cds.fasta

#3. adjustfasta
for f in *.fasta;
do 
    taxon=`echo $f | sed s/.pep.fasta//`;
    orthomclAdjustFasta $taxon $f 1;
done;

#4. orthomclFilterFasta

#5. makeblastdb -in goodProteins.fasta -dbtype 'prot'

#6. qsub gam.blastp.pbs

#7. cleanup database