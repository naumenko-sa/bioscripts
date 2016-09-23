#!/bin/bash

#exports gemini database to txt file to import to xls
#database schema: https://gemini.readthedocs.io/en/latest/content/database_schema.html#the-variants-table
#reports PASS only variants to the database

bname=`echo $1 | sed s/.db//`;
sample=$2

head -n1 $bname.tmp | awk '{for (i=1;i<=NF;i++) printf $i"\t";printf "exac_pLi\texac_mis_z\tomim\n"}' > $bname.txt
cat $bname.tmp | sed 1d > $bname.body;
cheo.add_exac_scores.pl $bname.body > $bname.body1;
cheo.add_omim_string.pl $bname.body1 >>  $bname.txt
#rm $bname.tmp $bname.body $bname.body1