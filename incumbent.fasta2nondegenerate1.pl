#!/usr/bin/perl

#get nondegenerate sites 
#prints 0-based site number (aa number),1,2,3 codon position
#requires 1 gapless triplets left and right
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
#dont look at first 1 and last 1 codons
for ($i=1;$i<length($aln[0])/3-1;$i++)
{
    $nondegenerate=0;
    @triplets;
    for ($j=0;$j<@aln;$j++)
    {
	$triplet = substr($aln[$j],$i*3,3);
	$nondegenerate += is_nondegenerate($triplet);
	push(@triplets,$triplet);
    }
    $dna='';
    for($j=0;$j<@aln;$j++)
    {
	$dna .= substr($aln[$j],$i*3-3,3*3); #left 1 triplets and right 1 triplets should not have gaps
    } 
    #print $dna."\n";
    if ($nondegenerate == @aln and $dna !~ /-/)
    {
	$sites++;
	$slice1 = '';
	$slice2 = '';
	$slice3 = '';
	for ($j=0;$j<@aln;$j++)
	{
	    $slice1 .= uc(substr($triplets[$j],0,1));
	    $slice2 .= uc(substr($triplets[$j],1,1));
	    $slice3 .= uc(substr($triplets[$j],2,1));
	}
	
	
	if( $slice1 =~ '([ATGC])\g1{4}' and $slice3  =~ '^([ATGC])\g1{4}' )
	{
	    print $i."\t".$slice1."\t".$slice2."\t".$slice3."\n";
	#    $syn++;
	}
    }
    while (@triplets>0) 
    {	
	shift @triplets;
    }
}

#print "Syn: ".$syn."\n";
#print "Sites: ".$sites."\n";

close(IN);

#64-9 = 55 triplets!
sub is_nondegenerate()
{
    $res = 1;
    if ($_[0] =~ /CT[ATGC]/i or $_[0] eq "TTA" or $_[0] eq "TTG" or $_[0] eq "TAA" or $_[0] eq "TAG" or $_[0] eq "TGA" or $_[0] =~ /-/)
    {
        $res=0;
    }
    return $res;
}