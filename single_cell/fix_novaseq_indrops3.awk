# adds N to the trucated reads and F to the quality values
# usage: gunzip -c R1.fq.gz | awk -v full_length=N fix_novaseq_indrops3.awk | bgzip > R1.corrected.fq.gz
# where full length = 8 for R2 and 14 for R4
{
    print $0;
    getline;
    l=length($0); 
    printf $0;
    for(i=0;i<full_length-l;i++)
	printf "N";
    printf "\n";
    getline;
    print $0;
    getline;
    printf $0;
    for(i=0;i<full_length-l;i++)
	printf "F";
    printf "\n";
}
