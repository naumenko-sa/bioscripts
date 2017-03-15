#!/bin/bash

ev=1e-5

usage()
{

    if [ $# -lt 3 ];
    then
	echo "Usage: gam_filter.sh left.fq right.fq base_of_contaminants.fa";
	exit;
    fi;

    bname=`echo $1 | sed s/.r1.fq//`

    nreads.sh  $1;
}

####################################################################################################
step1()
{

    echo "Step 1: trimming read names" >> gam_filter.log;
    date >> gam_filter.log;
    
    cat $1 | awk '{if(NR%4==1) print $1; else print $0;}' > $1.1
    cat $2 | awk '{if(NR%4==1) print $1; else print $0;}' > $2.1

    rm $1 $2

    mv $1.1 $1
    mv $2.1 $2
}

####################################################################################################

step2()
{

    echo `date` "Step 2: remove duplicates and low complexity reads" >> gam_filter.log;

    prinseq-lite.pl -fastq $1 -fastq2 $2 -log prinseq.log -derep 12345 -lc_method dust -lc_threshold 7 -out_good ${bname}.filt.fq

    mv ${bname}.filt.fq_1.fastq ${bname}.filt.r1.fq
    mv ${bname}.filt.fq_2.fastq ${bname}.filt.r2.fq

    cat *singletons.fastq > ${bname}.filt.single.fq

    rm *prinseq_bad* *singletons.fastq
    
    goodpairs=`cat prinseq.log | grep "Good sequences" | grep pairs | awk '{print $7}'`
    goodsingles=`cat prinseq.log | grep "Good sequences"  | grep singletons | awk '{print $9}' | sed s/,//g | awk '{sum+=$0}END{print sum}'`
    
    echo "Prinseq: "$goodpairs" goodpairs "$goodsingles" goodsingles" >> gam_filter.log;
}

####################################################################################################

step3()
{
    echo `date` "Step 3: find contamination" >> gam_filter.log;

    for f in {r1,r2,single};
    do
	fq2fa.sh ${bname}.filt.$f.fq > ${bname}.filt.$f.fa
    done

    makeblastdb -in $3 -dbtype 'nucl'

    for f in {r1,r2,single}
    do
	goblastn.sh ${bname}.filt.$f.fa $3 1 $ev
    done
    
    rm *.fasta.qual
}


####################################################################################################
step4()
{

    echo `date` "Step 4: remove reads" >> gam_filter.log;

    for f in {r1,r2,single}
    do
	cat ${bname}.filt.$f.fa_vs_$3.blastn.$ev | awk '{print $1}' >> ${bname}.cont.list.t
    done

    cat ${bname}.cont.list.t | sort | uniq > ${bname}.cont.list

    rm ${bname}.cont.list.t

    for f in {r1,r2,single}
    do
	get_readss_fastq.pl ${bname}.filt.$f.fq ${bname}.cont.list 0 > ${bname}.filt1.$f.fq
    done

    rm ${bname}.filt.*.fq ${bname}.filt.*.fa $3.* $bname.cont.list
    #*.blastn.${ev}

    for f in {r1,r2,single}
    do
	mv ${bname}.filt1.$f.fq ${bname}.filt.$f.fq
    done
    nreads.sh ${bname}.filt.r1.fq;
    nreads.sh ${bname}.filt.single.fq;
    
    echo "Blastn: "`cat ${bname}.filt.r1.fq.nreads`" goodpairs "`cat ${bname}.filt.single.fq.nreads`" goodsingles" >> gam_filter.log; 
}


####################################################################################################
step5()
{
    echo `date` "Step 5: merge overlapping reads" >> gam_filter.log;

    gofastqjoin.sh $bname.filt.r1.fq $bname.filt.r2.fq 20 $bname.merged.fq

    cat $bname.filt.single.fq >> $bname.merged.fqjoin
    rm $bname.filt.single.fq
    mv $bname.merged.fqjoin $bname.filt.single.fq

    rm $bname.filt.r?.fq

    mv $bname.merged.fqun1 $bname.filt.r1.fq
    mv $bname.merged.fqun2 $bname.filt.r2.fq

    for f in {r1,single}
    do
	nreads.sh ${bname}.filt.$f.fq;
    done
    
    echo "Fastqjoin: "`cat ${bname}.filt.r1.fq.nreads`" pairs "`cat ${bname}.filt.single.fq.nreads`" singles" >> gam_filter.log;
}
####################################################################################################
## step 5 fastqc and trimmomatic


usage "$@"
for f in {2,3,4,5}
do
    step$f "$@"
done

echo `date` end  >> gam_filter.log;