#!/usr/bin/perl

open(IN,$ARGV[0]);

while (<IN>)
{
    chomp;
    if ($_ !~ /CHROM/)
    {
	@ar = split("\t",$_);
	$id = $ar[0]."-".$ar[1]."-".$ar[2]."-".$ar[3];
	@ar2 = split(",",@ar[4]);
	print $id."\n";
	foreach $impact (@ar2)
	{
	    #print $impact."\n";
	    @ar3= split(/\|/,$impact);
	    if ($ar3[17] == '')
	    {
		push(@ar3,'NA');
	    }
	    print "exon",$ar3[9].":".$ar3[16].":".$ar3[17]."\n";
	}
    }
}

close(IN);