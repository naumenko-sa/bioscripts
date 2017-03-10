#!/usr/bin/perl

# parses vcf with VEP annotation - CSQ field

open(IN,$ARGV[0]);

while (<IN>)
{
    chomp;
    $impacts_string='';
    if ($_ !~ /CHROM/)
    {
	@ar = split("\t",$_);
	$id = $ar[0]."-".$ar[1]."-".$ar[2]."-".$ar[3];
	@ar2 = split(",",@ar[4]);
	#print $id."\n";
	foreach $impact (@ar2)
	{
	    #print $impact."\n";
	    @ar3= split(/\|/,$impact);
	    #if there is hvgsc
	    if ($ar3[16] ne '')
	    {
		if ($ar3[17] eq '')
		{
		    $ar3[17] = 'NA';
		}
		$impacts_string .= "exon".$ar3[9].":".$ar3[16].":".$ar3[17].",";
	    }
	}
	print $id."\t".substr($impacts_string,0,-1)."\n";
    }
}

close(IN);