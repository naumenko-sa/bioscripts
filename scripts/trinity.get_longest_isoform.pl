#!/usr/bin/perl

open(IN,$ARGV[0]);

%locus_len;
%locus_dna;

while(<IN>)
{
    $name = $_;
    $dna = <IN>;
    
    @ar = split(" ",$name);
    @ar1 = split ("_",$ar[0]);
    $locus_name = $ar1[0]."_".$ar1[1];
    
    @ar1 = split("=",$ar[1]);
    $len = $ar1[1];
    if (exists($locus_len{$locus_name}))
    {
	if ($len > $locus_len{$locus_name})
	{
	    $locus_len{$locus_name} = $len;
	    $locus_dna{$locus_name} = $dna;
	}
    }
    else
    {
	$locus_len{$locus_name}=$len;
	$locus_dna{$locus_name}=$dna;
    }
    
}

foreach $key (sort keys(%locus_len))
{
    print $key."_len=".$locus_len{$key}."\n";
    print $locus_dna{$key};
}

close(IN);