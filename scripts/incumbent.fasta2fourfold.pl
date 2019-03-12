#!/usr/bin/perl

#get 4-fold synonymous sites 
#prints 0-based site number (aa number), 3 codon position slice
#requires 9 gapless triplets left and right
#requires consevative left (+2, this site) and right(+4, next site) slice

@aln;

open(IN,$ARGV[0]);
while(<IN>)
{
    chomp;
    $name=$_;
    $dna = <IN>;
    chomp $dna;
    push (@aln,$dna);
}

$par=0;
$syn=0;
$sites=0;
#dont look at first 9 and last 9 codons
for ($i=9;$i<length($aln[0])/3-9;$i++)
{
    $fourfold=0;
    @triplets;
    for ($j=0;$j<@aln;$j++)
    {
	$triplet = substr($aln[$j],$i*3,4);
	$fourfold += is_4fold($triplet);
	push(@triplets,$triplet);
    }
    $dna='';
    for($j=0;$j<@aln;$j++)
    {
	$dna .= substr($aln[$j],$i*3-27,19*3); #left 9 triplets and right 9 triplets should not have gaps
    } 
    #print $dna."\n";
    if ($fourfold == @aln and $dna !~ /-/)
    {
	$sites++;
	$slice1 = '';
	$slice2 = '';
	$slice3 = '';
	$slice4 = '';
	for ($j=0;$j<@aln;$j++)
	{
	    #$slice1 .= uc(substr($triplets[$j],0,1));
	    $slice2 .= uc(substr($triplets[$j],1,1));
	    $slice3 .= uc(substr($triplets[$j],2,1));
	    $slice4 .= uc(substr($triplets[$j],3,1));
	}
	
	#depends on the number of species
	if ($slice2 =~ '([ATGC])\g1{4}' and $slice4 =~ '^([ATGC])\g1{4}')
	{
	    #print $i."\t".$slice1."\t".$slice2."\t".$slice3."\t".$slice4."\n";
	    print $i."\t".$slice3."\n";
	}
    }
    while (@triplets>0) 
    {	
	shift @triplets;
    }
}

close(IN);

#8 codons = 32 triplets!
sub is_4fold()
{
    $res=0;
    if ( $_[0] =~ /TC[ATGC]{2}/i or 
	$_[0] =~ /CT[ATGC]{2}/i or 
	$_[0] =~ /CC[ATGC]{2}/i or 
	$_[0] =~ /CG[ATGC]{2}/i or 
        $_[0] =~ /AC[ATGC]{2}/i or 
        $_[0] =~ /GT[ATGC]{2}/i or 
        $_[0] =~ /GC[ATGC]{2}/i or 
        $_[0] =~ /GG[ATGC]{2}/i)
    {
	$res=1;
	#print $_[0]."\n";
    }
    return $res;
}