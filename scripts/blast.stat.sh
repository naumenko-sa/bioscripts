#/bin/bash

#takes one longest hit

echo "Unique_hits:" `cat $1 | awk '{print $1}' | sort |uniq | wc -l`
cat $1 | awk '{if ($3>=200) print $0;}' | sort -k 1,3 -nr > $1.200
cat $1.200 | awk '{if ( ($1 == prev1) &&  ($2 == prev2)) next;else {prev1 = $1; prev2=$2;print $0;}}' > $1.uniq
echo ">200 hits:" `cat $1.uniq  | wc -l`
echo "total len:" `cat $1.uniq | awk '{sum+=$3} END{print sum}'`
echo "avr ID:" `cat $1.uniq  |  awk '{id+=$4} END{print id/NR}'`

#gam=gam4.2;for f in `ls *_vs_$gam.fasta*.1e-10`;do echo $f >> $gam.stat ;goblastn_stat.sh $f >>$gam.stat;done;
#gam=gam4.1;cat $gam.stat | awk '{if (NR %5 == 1 || NR%5 ==0) print $0;}' | sed s/.fasta_vs_$gam.fasta.blastn.1e-10// | sed s/avr\ ID:// | awk '{first=$0;getline;print first"\t"$0;}'> $gam.dist
#gam=gam4.1;for f in `cat species.order`;do egrep "^${f}\s" $gam.dist;done;