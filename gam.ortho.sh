#!/bin/bash

# run orthomcl to build cogs for a group of species in species
# manual: http://orthomcl.org/common/downloads/software/v2.0/UserGuide.txt
# this is not a script, it is a history or log file

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
for taxon in `cat specie_s`;
do 
    orthomclAdjustFasta $taxon $taxon.cd-hit.pep.fasta 1;
done;
#do this for cds as well, because they are necessary for reverse translation

#4. orthomclFilterFasta

#5. makeblastdb -in goodProteins.fasta -dbtype 'prot'

#6. qsub gam.blastp.pbs

#7. cleanup database

#8. blast parser
orthomclBlastParser goodProteins.fasta_vs_goodProteins.fasta.blastp.1e-5 compliantFasta/ > similarSequences.txt

#9. load similar sequences
orthomclLoadBlast ../config_file similarSequences.txt

#10. orthomcl pairs
orthomclPairs ../config_file log_file cleanup=yes

#11. dump pair files
orthomclDumpPairsFiles ../config_file

#12 mcl
mcl mclInput --abc -I 1.5 -o mclOutput

#13
orthomclMclToGroups cog 00001 < mclOutput > groups.txt

#14 - extract 1:1 orthologs for 6 species
gam.orthomcl.make_alignments.pl groups.txt goodProteins.fasta 6

#15 align protein sequences
for f in *.fasta;do alignment.muscle.sh $f;done;

#16 select conservative alignments
for f in *.align.fasta;do gaps=`alignment.count_gaps.sh $f`;if [ $gaps -le 10 ];then mv $f conservative/;fi;done;

#17. reverse translate alignments
alignment.pep2dna.sh