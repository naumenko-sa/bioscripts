$!/bin/bash

#cat knownCanonical.exonNuc.fa.noempty | awk '{name=$0;getline;dna=$0; if ((name ~ /mm9/) || (name ~ /rn4/) || (name ~ /oryCun2/) || (name ~ /ochPri2/ )) {print name;print dna;}}' > dist6.fa
#cat $1.fa  | awk 'BEGIN{dna="";name=""}{if (NR%8==1) {name=name$0;getline;dna=dna$0;}else{getline;}}END{print ">mm9";print dna}' > $1.fasta

echo "synonymous divergence"

fasta2fourfold.pl $1.fasta > $1.4fold

wc -l $1.4fold 

cat $1.4fold | awk '{print $3}' | awk -F "" '{if ($1 != $2) print $0}' | wc -l

cat $1.4fold | awk '{print $3}' | awk -F "" '{if ($1 != $2 && $1 == $3 && $2 == $4) print $0}' | wc -l
cat $1.4fold | awk '{print $3}' | awk -F "" '{if ($1 != $2 && $1 == $4 && $2 == $3) print $0}' | wc -l

cat $1.4fold | awk '{print $3}' | awk -F "" '{if ($1==$2 && $3!=$4 && $1==$3) print $0;}' | wc -l
cat $1.4fold | awk '{print $3}' | awk -F "" '{if ($1==$2 && $3!=$4 && $1==$4) print $0;}' | wc -l
                                      


echo "nonsynonymous divergence"

fasta2nondegenerate.pl $1.fasta > $1.nondeg
cat $1.nondeg | awk '{if (($1=="AAAA" || $1=="CCCC" || $1=="GGGG" || $1=="TTTT") && ($2 != "AAAA" && $2 != "TTTT" && $2 != "CCCC" && $2 != "GGGG")) print $0;}' > $1.nondeg1
cat $1.nondeg | awk '{if (($2=="AAAA" || $2=="CCCC" || $2=="GGGG" || $2=="TTTT") && ($1 != "AAAA" && $1 != "TTTT" && $1 != "CCCC" && $1 != "GGGG")) print $0;}' > $1.nondeg2

cat $1.nondeg1 | awk '{print $2}' | awk -F "" '{if ($1!=$2) print $0}' | wc -l
cat $1.nondeg2 | awk '{print $1}' | awk -F "" '{if ($1!=$2) print $0}' | wc -l 

wc -l $1.nondeg1
wc -l $1.nondeg2

echo "parallel nonsynonymous divergence"

cat $1.nondeg1 | awk '{print $2}' | awk -F "" '{if ($1!=$2 && $3!=$4 && $1==$3) print $0}' | wc -l
cat $1.nondeg1 | awk '{print $2}' | awk -F "" '{if ($1!=$2 && $3!=$4 && $1==$4) print $0}' | wc -l

cat $1.nondeg2 | awk '{print $1}' | awk -F "" '{if ($1!=$2 && $3!=$4 && $1==$4) print $0}' | wc -l
cat $1.nondeg2 | awk '{print $1}' | awk -F "" '{if ($1!=$2 && $3!=$4 && $1==$3) print $0}' | wc -l

cat $1.nondeg2 | awk '{print $1}' | awk -F "" '{if ($1==$2 && $3!=$4 && $1==$3) print $0}' | wc -l
cat $1.nondeg2 | awk '{print $1}' | awk -F "" '{if ($1==$2 && $3!=$4 && $1==$4) print $0}' | wc -l
cat $1.nondeg1 | awk '{print $2}' | awk -F "" '{if ($1==$2 && $3!=$4 && $1==$4) print $0}' | wc -l
cat $1.nondeg1 | awk '{print $2}' | awk -F "" '{if ($1==$2 && $3!=$4 && $1==$3) print $0}' | wc -l